"""Instrumented copy of the reference Newton-PCG step, with optional inexact Newton."""

from __future__ import annotations

import time

import numpy as np
from scipy import fft
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg


def advance_profiled(sim, dt_s, source_W_m3=0.0, boundary=None, linear_policy="fixed",
                     fft_workers=1, preconditioner_method="fft"):
    """Advance one step and return a detailed profile without changing the model.

    ``fixed`` uses the reference 1e-8 (or configured) PCG tolerance. ``ew`` uses
    a conservative Eisenstat-Walker forcing term for early Newton iterations;
    the nonlinear acceptance test and final energy accounting are unchanged.
    """
    from implicit_enthalpy_fvm import BoundarySet

    if dt_s <= 0.0 or not np.isfinite(dt_s):
        raise ValueError("dt_s must be finite and positive")
    if linear_policy not in {"fixed", "ew", "ew_safe"}:
        raise ValueError(f"unknown linear tolerance policy: {linear_policy}")
    if preconditioner_method not in {"fft", "amg"}:
        raise ValueError(f"unknown preconditioner: {preconditioner_method}")
    boundary = boundary or BoundarySet()
    boundary.validate(sim.ndim)
    source = np.broadcast_to(np.asarray(source_W_m3, dtype=np.float64), sim.shape)
    if np.any(~np.isfinite(source)):
        raise ValueError("source must be finite")

    totals = {
        "property_update_seconds": 0.0, "property_update_calls": 0,
        "face_resistance_seconds": 0.0, "face_resistance_calls": 0,
        "face_flux_seconds": 0.0, "face_flux_calls": 0,
        "matrix_vector_seconds": 0.0, "matrix_vector_calls": 0,
        "fft_forward_seconds": 0.0, "fft_forward_calls": 0,
        "fft_inverse_seconds": 0.0, "fft_inverse_calls": 0,
        "pcg_seconds": 0.0, "pcg_calls": 0,
    }
    originals = []

    def time_method(owner, name, seconds_key, calls_key):
        original = getattr(owner, name)

        def wrapped(*args, **kwargs):
            started = time.perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                totals[seconds_key] += time.perf_counter() - started
                totals[calls_key] += 1

        setattr(owner, name, wrapped)
        originals.append((owner, name, original))

    # Inclusive high-level times are paired with lower-level timings; do not sum
    # these columns together because residual evaluation contains those calls.
    for name in ("enthalpy_density_J_m3", "volumetric_capacity_J_m3_K", "conductivity_W_m_K"):
        time_method(sim.materials, name, "property_update_seconds", "property_update_calls")
    time_method(sim, "_face_conductances_density", "face_resistance_seconds", "face_resistance_calls")
    time_method(sim, "_net_conduction_density", "face_flux_seconds", "face_flux_calls")
    time_method(sim, "_apply_laplacian", "matrix_vector_seconds", "matrix_vector_calls")

    rfft_original, irfft_original = fft.rfftn, fft.irfftn

    def timed_rfft(*args, **kwargs):
        started = time.perf_counter()
        try:
            return rfft_original(*args, **kwargs)
        finally:
            totals["fft_forward_seconds"] += time.perf_counter() - started
            totals["fft_forward_calls"] += 1

    def timed_irfft(*args, **kwargs):
        started = time.perf_counter()
        try:
            return irfft_original(*args, **kwargs)
        finally:
            totals["fft_inverse_seconds"] += time.perf_counter() - started
            totals["fft_inverse_calls"] += 1

    fft.rfftn, fft.irfftn = timed_rfft, timed_irfft
    try:
        records, iterations = _advance_profiled_core(
            sim, dt_s, source, boundary, linear_policy, totals, fft_workers,
            preconditioner_method,
        )
    finally:
        fft.rfftn, fft.irfftn = rfft_original, irfft_original
        for owner, name, original in reversed(originals):
            setattr(owner, name, original)
    return records, {"components": totals, "newton_iterations": iterations}


def _advance_profiled_core(sim, dt_s, source, boundary, linear_policy, totals, fft_workers,
                           preconditioner_method):
    h_old = sim.enthalpy_density_J_m3
    t = sim.temperature_K.copy()
    total_cg_iterations = 0
    max_cg_relative_residual = 0.0
    nonlinear_residual = np.inf
    previous_residual = None
    iteration_log = []
    step_start = time.perf_counter()

    for nonlinear_iteration in range(1, sim.max_newton + 1):
        iter_start = time.perf_counter()
        begin = totals.copy()
        residual, gfaces, gbfaces, fields = sim._residual(t, h_old, dt_s, source, boundary)
        nonlinear_residual = sim._relative_residual(
            residual, h_old, dt_s * source, dt_s * fields["conduction"]
        )
        if nonlinear_residual <= sim.nonlinear_rtol:
            iteration_log.append({
                "newton_iteration": nonlinear_iteration,
                "relative_residual": nonlinear_residual,
                "pcg_iterations": 0,
                "pcg_relative_residual": 0.0,
                "linear_rtol": 0.0,
                "line_search_trials": 0,
                "accepted_damping": 0.0,
                "wall_seconds": time.perf_counter() - iter_start,
                "already_converged": True,
            })
            break

        capacity = sim.materials.volumetric_capacity_J_m3_K(t, sim.ag_fraction)

        def matvec(flat):
            delta_t = flat.reshape(sim.shape)
            out = capacity * delta_t + dt_s * sim._apply_laplacian(delta_t, gfaces, gbfaces)
            return out.ravel()

        operator = LinearOperator((sim.ag_fraction.size, sim.ag_fraction.size), matvec=matvec, dtype=np.float64)
        setup_start = time.perf_counter()
        mean_capacity = float(np.mean(capacity))
        spectrum_shape = tuple(sim.shape[:-1]) + (sim.shape[-1] // 2 + 1,)
        spectrum = np.full(spectrum_shape, mean_capacity, dtype=np.float64)
        for axis, g in enumerate(gfaces):
            n_axis = sim.shape[axis]
            if axis == sim.ndim - 1:
                frequencies = np.arange(n_axis // 2 + 1, dtype=np.float64) / n_axis
            else:
                frequencies = np.arange(n_axis, dtype=np.float64) / n_axis
            shape = [1] * sim.ndim
            shape[axis] = frequencies.size
            lam = 2.0 * float(np.mean(g)) * (1.0 - np.cos(2.0 * np.pi * frequencies))
            spectrum += dt_s * lam.reshape(shape)
        preconditioner_setup_seconds = time.perf_counter() - setup_start

        if linear_policy == "fixed":
            linear_rtol = sim.linear_rtol
        elif previous_residual is None and linear_policy == "ew":
            linear_rtol = min(0.25, max(sim.linear_rtol, 0.05))
        elif previous_residual is None:
            linear_rtol = min(0.05, max(sim.linear_rtol, 0.005))
        elif linear_policy == "ew_safe":
            ratio = max(nonlinear_residual / max(previous_residual, 1.0e-30), 1.0e-16)
            linear_rtol = min(0.1, max(sim.linear_rtol, 0.005, 0.5 * ratio ** 1.5))
        else:
            ratio = max(nonlinear_residual / max(previous_residual, 1.0e-30), 1.0e-16)
            linear_rtol = min(0.5, max(sim.linear_rtol, 0.02, 0.9 * ratio ** 1.5))

        amg_setup_seconds = 0.0
        amg_levels = None
        if preconditioner_method == "fft":
            preconditioner = LinearOperator(
                (sim.ag_fraction.size, sim.ag_fraction.size),
                matvec=lambda flat: fft.irfftn(
                    fft.rfftn(flat.reshape(sim.shape), workers=fft_workers) / spectrum,
                    s=sim.shape, workers=fft_workers,
                ).ravel(), dtype=np.float64,
            )
        else:
            # Pilot-only sparse Galerkin hierarchy; build the exact symmetric
            # seven-point Jacobian from the same face conductances as the FVM.
            import pyamg

            amg_started = time.perf_counter()
            n_cells = sim.ag_fraction.size
            index = np.arange(n_cells, dtype=np.int64).reshape(sim.shape)
            diagonal = capacity.ravel().copy()
            rows = []
            columns = []
            values = []
            for axis, face_g in enumerate(gfaces):
                diagonal += dt_s * (face_g + np.roll(face_g, 1, axis=axis)).ravel()
                neighbor = np.roll(index, -1, axis=axis)
                edge_values = -dt_s * face_g
                rows.append(index.ravel())
                columns.append(neighbor.ravel())
                values.append(edge_values.ravel())
                rows.append(neighbor.ravel())
                columns.append(index.ravel())
                values.append(edge_values.ravel())
            for gb in gbfaces.values():
                diagonal += dt_s * gb.ravel()
            rows.append(index.ravel())
            columns.append(index.ravel())
            values.append(diagonal)
            jacobian = sparse.coo_matrix(
                (np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
                shape=(n_cells, n_cells),
            ).tocsr()
            jacobian.eliminate_zeros()
            hierarchy = pyamg.smoothed_aggregation_solver(jacobian, max_levels=10, max_coarse=25)
            preconditioner = hierarchy.aspreconditioner(cycle="V")
            amg_levels = len(hierarchy.levels)
            amg_setup_seconds = time.perf_counter() - amg_started
        cg_count = [0]

        def cg_callback(_x):
            cg_count[0] += 1

        cg_start = time.perf_counter()
        delta, info = cg(
            operator, -residual.ravel(), M=preconditioner, rtol=linear_rtol,
            atol=0.0, maxiter=sim.max_cg, callback=cg_callback,
        )
        cg_seconds = time.perf_counter() - cg_start
        totals["pcg_seconds"] += cg_seconds
        totals["pcg_calls"] += 1
        total_cg_iterations += cg_count[0]
        linear_residual = float(np.linalg.norm(matvec(delta) + residual.ravel())
                                / max(np.linalg.norm(residual), 1.0e-30))
        max_cg_relative_residual = max(max_cg_relative_residual, linear_residual)
        if info != 0:
            raise RuntimeError(
                f"PCG did not converge (info={info}, iterations={cg_count[0]}, rel={linear_residual:.3e})"
            )

        accepted = False
        damping = 1.0
        line_search_trials = 0
        line_search_start = time.perf_counter()
        for line_search_trials in range(1, 33):
            candidate = t + damping * delta.reshape(sim.shape)
            candidate_max = sim.materials.ti_Tmax_K if np.any(sim.ag_fraction < 1.0 - 1.0e-14) else 10000.0
            if float(candidate.min()) < 100.0 or float(candidate.max()) > candidate_max + 1.0e-8:
                damping *= 0.5
                continue
            new_residual, _, _, new_fields = sim._residual(candidate, h_old, dt_s, source, boundary)
            candidate_norm = sim._relative_residual(
                new_residual, h_old, dt_s * source, dt_s * new_fields["conduction"]
            )
            if candidate_norm < nonlinear_residual or candidate_norm <= sim.nonlinear_rtol:
                previous_residual = nonlinear_residual
                t = candidate
                nonlinear_residual = candidate_norm
                accepted = True
                break
            damping *= 0.5
        line_search_seconds = time.perf_counter() - line_search_start
        if not accepted:
            raise RuntimeError(
                f"Newton line search failed at iteration {nonlinear_iteration}; "
                f"relative residual={nonlinear_residual:.3e}"
            )

        component_delta = {key: totals[key] - begin[key] for key in begin}
        iteration_log.append({
            "newton_iteration": nonlinear_iteration,
            "relative_residual": nonlinear_residual,
            "pcg_iterations": cg_count[0],
            "pcg_relative_residual": linear_residual,
            "linear_rtol": linear_rtol,
            "line_search_trials": line_search_trials,
            "accepted_damping": damping,
            "preconditioner": preconditioner_method,
            "amg_setup_seconds": amg_setup_seconds,
            "amg_levels": amg_levels,
            "preconditioner_setup_seconds": preconditioner_setup_seconds,
            "line_search_seconds": line_search_seconds,
            "wall_seconds": time.perf_counter() - iter_start,
            "already_converged": False,
            **component_delta,
        })
        if nonlinear_residual <= sim.nonlinear_rtol:
            break
    else:
        raise RuntimeError(
            f"Newton did not converge in {sim.max_newton} iterations; relative residual={nonlinear_residual:.3e}"
        )

    final_residual, _, _, final_fields = sim._residual(t, h_old, dt_s, source, boundary)
    nonlinear_residual = sim._relative_residual(
        final_residual, h_old, dt_s * source, dt_s * final_fields["conduction"]
    )
    if nonlinear_residual > max(sim.nonlinear_rtol * 1.01, 1.0e-12):
        raise RuntimeError(f"post-step nonlinear residual is too high: {nonlinear_residual:.3e}")
    liq = sim.materials.liquid_fraction(t)
    active = sim.ag_fraction > 1.0e-14
    if np.any(liq[active] < -1.0e-12) or np.any(liq[active] > 1.0 + 1.0e-12):
        raise RuntimeError("nonphysical Ag liquid fraction; no clipping is applied")
    if np.any(sim.ag_fraction < 1.0 - 1.0e-14) and float(t.max()) > sim.materials.ti_Tmax_K + 1.0e-8:
        raise RuntimeError("Ti3SiC2 temperature exceeded the audited candidate range")

    source_rate_W = float(np.sum(source) * sim.volume_m3)
    out_rate_W = sim._boundary_out_rate_W(t, final_fields["conductivity"], boundary)
    source_increment = dt_s * source_rate_W
    out_increment = dt_s * out_rate_W
    storage_change = float(np.sum(final_fields["enthalpy"] - h_old) * sim.volume_m3)
    sim.temperature_K = t
    sim.enthalpy_density_J_m3 = final_fields["enthalpy"]
    sim.time_s += dt_s
    sim.source_energy_J += source_increment
    sim.boundary_energy_out_J += out_increment
    sim.step_count += 1
    storage_total = float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3)
    energy_residual = (storage_total - sim.initial_stored_enthalpy_J
                       - sim.source_energy_J + sim.boundary_energy_out_J)
    active_liquid = liq[active]
    record = {
        "time_s": sim.time_s, "dt_s": dt_s, "nonlinear_iterations": len(iteration_log),
        "pcg_iterations_total": total_cg_iterations,
        "pcg_relative_residual_max": max_cg_relative_residual,
        "nonlinear_relative_residual": nonlinear_residual,
        "max_temperature_K": float(t.max()),
        "mean_ag_temperature_K": float(np.sum(sim.ag_fraction * t) / max(np.sum(sim.ag_fraction), 1.0e-30)),
        "max_ag_liquid_fraction": float(np.max(active_liquid)) if active_liquid.size else 0.0,
        "ag_liquid_volume_fraction": float(np.sum(sim.ag_fraction * liq) / max(np.sum(sim.ag_fraction), 1.0e-30)),
        "source_energy_step_J": source_increment, "boundary_energy_out_step_J": out_increment,
        "storage_change_step_J": storage_change, "total_stored_enthalpy_J": storage_total,
        "cumulative_source_energy_J": sim.source_energy_J,
        "cumulative_boundary_out_J": sim.boundary_energy_out_J,
        "energy_residual_J": energy_residual, "step_wall_seconds": time.perf_counter() - step_start,
    }
    return record, iteration_log
