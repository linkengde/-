#!/usr/bin/env python3
"""Analysis routines for the C-frozen Ag64 rho100 follow-up plan.

This module performs no dynamics.  Its metric choices are versioned in the
execution config so C can review them before issuing a run release.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Iterable, Sequence

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from scipy.special import sph_harm_y


Q6_NEIGHBOR_CUTOFF_A = 3.8
CNA_NEIGHBOR_CUTOFF_A = 3.8
RDF_BIN_WIDTH_A = 0.05
WINDOW_STEPS = 4000
THERMAL_LOG_INTERVAL_STEPS = 20
FRAME_INTERVAL_STEPS = 100


class CriteriaError(RuntimeError):
    pass


@dataclass(frozen=True)
class FrameMetrics:
    step: int
    temperature_K: float
    potential_energy_eV_per_atom: float
    global_q6: float
    crystalline_fraction: float
    largest_crystalline_cluster_atoms: int
    coordination_number: float
    rdf_first_peak_A: float | None
    rdf_first_minimum_A: float | None
    minimum_pair_distance_A: float
    rdf_r_A: tuple[float, ...]
    rdf_g: tuple[float, ...]
    stress_eV_A3: tuple[float, ...]


def _finite_array(value: object, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise CriteriaError(f"{name} contains non-finite values")
    return array


def _periodic_neighbor_data(atoms: Atoms, cutoff_A: float):
    i, j, shifts, distances = neighbor_list("ijSd", atoms, cutoff_A, self_interaction=False)
    cell = np.asarray(atoms.cell, dtype=float)
    vectors = np.asarray(atoms.positions)[j] + np.asarray(shifts) @ cell - np.asarray(atoms.positions)[i]
    return i.astype(int), j.astype(int), shifts.astype(int), np.asarray(distances, dtype=float), vectors


def global_q6(atoms: Atoms, cutoff_A: float = Q6_NEIGHBOR_CUTOFF_A) -> float:
    """Global Steinhardt Q6, retaining periodic image offsets from ASE."""
    i, _j, _s, distances, vectors = _periodic_neighbor_data(atoms, cutoff_A)
    n_atoms = len(atoms)
    if n_atoms == 0:
        raise CriteriaError("Q6 requires a nonempty system")
    qlm = np.zeros((n_atoms, 13), dtype=np.complex128)
    degree = np.zeros(n_atoms, dtype=int)
    l = 6
    theta = np.arccos(np.clip(vectors[:, 2] / distances, -1.0, 1.0))
    phi = np.arctan2(vectors[:, 1], vectors[:, 0])
    for column, m in enumerate(range(-l, l + 1)):
        qlm[:, column] = np.bincount(
            i,
            weights=np.real(sph_harm_y(l, m, theta, phi)),
            minlength=n_atoms,
        ) + 1j * np.bincount(
            i,
            weights=np.imag(sph_harm_y(l, m, theta, phi)),
            minlength=n_atoms,
        )
    degree[:] = np.bincount(i, minlength=n_atoms)
    if np.any(degree == 0):
        raise CriteriaError("Q6 found an atom with no neighbors inside the locked cutoff")
    qlm /= degree[:, None]
    mean_qlm = qlm.mean(axis=0)
    result = math.sqrt(4.0 * math.pi / 13.0 * float(np.sum(np.abs(mean_qlm) ** 2)))
    if not math.isfinite(result):
        raise CriteriaError("Q6 is non-finite")
    return result


def _longest_path(adjacency: dict[tuple[int, int, int, int], set[tuple[int, int, int, int]]]) -> int:
    """Longest simple path in the tiny common-neighbor graph (at most four sites)."""
    best = 0

    def visit(node, seen):
        nonlocal best
        best = max(best, len(seen))
        for nxt in adjacency[node]:
            if nxt not in seen:
                visit(nxt, seen | {nxt})

    for node in adjacency:
        visit(node, {node})
    return best


def cna_atom_classes(atoms: Atoms, cutoff_A: float = CNA_NEIGHBOR_CUTOFF_A) -> tuple[list[str], list[set[int]]]:
    """Honeycutt-Andersen common-neighbor analysis for FCC/HCP Ag.

    Neighbor identities include lattice image offsets.  FCC atoms have 12
    (4,2,1) bonds; HCP atoms have six (4,2,1) and six (4,2,2) bonds.
    This is a topology-based CNA implementation, not a PTM fit.
    """
    i, j, shifts, _distances, _vectors = _periodic_neighbor_data(atoms, cutoff_A)
    n = len(atoms)
    by_center: list[list[tuple[int, tuple[int, int, int]]]] = [[] for _ in range(n)]
    for a, b, shift in zip(i, j, shifts):
        by_center[int(a)].append((int(b), tuple(int(x) for x in shift)))

    zero = (0, 0, 0)
    classes: list[str] = []
    crystal_graph: list[set[int]] = [set() for _ in range(n)]
    for center in range(n):
        center_sites = {(atom, shift) for atom, shift in by_center[center]}
        signatures: list[tuple[int, int, int]] = []
        for neighbor_atom, neighbor_shift in by_center[center]:
            neighbor_site = (neighbor_atom, neighbor_shift)
            around_neighbor = {
                (other_atom, tuple(neighbor_shift[k] + other_shift[k] for k in range(3)))
                for other_atom, other_shift in by_center[neighbor_atom]
            }
            common = center_sites & around_neighbor
            common.discard((center, zero))
            common.discard(neighbor_site)
            induced: dict[tuple[int, int, int, int], set[tuple[int, int, int, int]]] = {}
            normalized_common = {(atom, *shift) for atom, shift in common}
            for common_atom, common_shift in common:
                node = (common_atom, *common_shift)
                adjacent = set()
                for other_atom, other_shift in by_center[common_atom]:
                    absolute = tuple(common_shift[k] + other_shift[k] for k in range(3))
                    candidate = (other_atom, *absolute)
                    if candidate in normalized_common:
                        adjacent.add(candidate)
                induced[node] = adjacent
            bond_count = sum(len(v) for v in induced.values()) // 2
            signatures.append((len(common), bond_count, _longest_path(induced) if induced else 0))
            if neighbor_atom != center:
                crystal_graph[center].add(neighbor_atom)

        count_421 = signatures.count((4, 2, 1))
        count_422 = signatures.count((4, 2, 2))
        if len(signatures) == 12 and count_421 == 12:
            classes.append("FCC")
        elif len(signatures) == 12 and count_421 == 6 and count_422 == 6:
            classes.append("HCP")
        else:
            classes.append("OTHER")

    crystalline = {idx for idx, name in enumerate(classes) if name in {"FCC", "HCP"}}
    seen: set[int] = set()
    largest = 0
    for start in crystalline:
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            for nxt in crystal_graph[node] & crystalline:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        largest = max(largest, size)
    return classes, crystal_graph


def cna_metrics(atoms: Atoms, cutoff_A: float = CNA_NEIGHBOR_CUTOFF_A) -> tuple[float, int, list[str]]:
    classes, _graph = cna_atom_classes(atoms, cutoff_A)
    fraction = sum(name in {"FCC", "HCP"} for name in classes) / len(classes)
    crystalline = {idx for idx, name in enumerate(classes) if name in {"FCC", "HCP"}}
    i, j, _s, _d, _v = _periodic_neighbor_data(atoms, cutoff_A)
    graph = [set() for _ in classes]
    for a, b in zip(i, j):
        if a != b and a in crystalline and b in crystalline:
            graph[int(a)].add(int(b))
    seen: set[int] = set()
    largest = 0
    for start in crystalline:
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            for nxt in graph[node]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        largest = max(largest, size)
    return fraction, largest, classes


def _half_shortest_cell_height(cell: np.ndarray) -> float:
    inv = np.linalg.inv(cell)
    heights = [1.0 / np.linalg.norm(inv[:, axis]) for axis in range(3)]
    return min(heights) / 2.0


def rdf_metrics(atoms: Atoms, bin_width_A: float = RDF_BIN_WIDTH_A) -> dict:
    """Finite-cell RDF through half the shortest cell height."""
    rmax = _half_shortest_cell_height(np.asarray(atoms.cell, dtype=float))
    edges = np.arange(0.0, rmax + bin_width_A, bin_width_A)
    if edges[-1] > rmax:
        edges[-1] = rmax
    if len(edges) < 3:
        raise CriteriaError("cell is too small for the RDF radial grid")
    i, j, _s, distances, _v = _periodic_neighbor_data(atoms, rmax)
    counts, _ = np.histogram(distances, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    shell = 4.0 * math.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    rho = len(atoms) / atoms.get_volume()
    g = counts / (len(atoms) * rho * shell)
    peak_mask = (centers >= 2.6) & (centers <= 3.2)
    peak_idx_candidates = np.flatnonzero(peak_mask)
    peak_idx = int(peak_idx_candidates[np.argmax(g[peak_mask])]) if len(peak_idx_candidates) else None
    valley_idx = None
    if peak_idx is not None:
        valley_mask = (centers >= 3.4) & (centers <= 4.2) & (np.arange(len(centers)) > peak_idx)
        valley_idx_candidates = np.flatnonzero(valley_mask)
        if len(valley_idx_candidates):
            valley_idx = int(valley_idx_candidates[np.argmin(g[valley_mask])])
    coordination = float(np.count_nonzero(distances <= centers[valley_idx]) / len(atoms)) if valley_idx is not None else math.nan
    min_distance = float(np.min(distances)) if len(distances) else math.inf
    return {
        "r_A": centers.tolist(),
        "g": g.tolist(),
        "first_peak_A": float(centers[peak_idx]) if peak_idx is not None else None,
        "first_minimum_A": float(centers[valley_idx]) if valley_idx is not None else None,
        "coordination_number": coordination,
        "minimum_pair_distance_A": min_distance,
        "neighbor_count": int(len(i)),
    }


def frame_metrics(atoms: Atoms, step: int, temperature_K: float, potential_energy_eV: float, stress_eV_A3: Sequence[float]) -> FrameMetrics:
    pe = float(potential_energy_eV) / len(atoms)
    q6 = global_q6(atoms)
    crystal_fraction, largest, _classes = cna_metrics(atoms)
    rdf = rdf_metrics(atoms)
    stress = _finite_array(stress_eV_A3, "stress").reshape(-1)
    return FrameMetrics(
        step=int(step),
        temperature_K=float(temperature_K),
        potential_energy_eV_per_atom=pe,
        global_q6=q6,
        crystalline_fraction=crystal_fraction,
        largest_crystalline_cluster_atoms=int(largest),
        coordination_number=float(rdf["coordination_number"]),
        rdf_first_peak_A=rdf["first_peak_A"],
        rdf_first_minimum_A=rdf["first_minimum_A"],
        minimum_pair_distance_A=float(rdf["minimum_pair_distance_A"]),
        rdf_r_A=tuple(rdf["r_A"]),
        rdf_g=tuple(rdf["g"]),
        stress_eV_A3=tuple(float(x) for x in stress),
    )


def integrated_autocorrelation_time(values: Sequence[float]) -> float:
    """Initial-positive-sequence estimate of tau_int, in sample intervals."""
    x = _finite_array(values, "autocorrelation series")
    n = len(x)
    if n < 4:
        return 0.5
    x = x - x.mean()
    var = float(np.dot(x, x) / n)
    if var <= np.finfo(float).eps:
        return 0.5
    tau = 0.5
    for lag in range(1, n):
        rho = float(np.dot(x[:-lag], x[lag:]) / ((n - lag) * var))
        if rho <= 0:
            break
        tau += rho
    return max(0.5, tau)


def block_bootstrap_mean_ci(values: Sequence[float], confidence: float = 0.95, replicates: int = 2000) -> tuple[float, float, int]:
    """Circular moving-block bootstrap; deterministic analysis-only seed."""
    x = _finite_array(values, "bootstrap series")
    if len(x) < 4:
        raise CriteriaError("at least four samples are needed for block bootstrap")
    block = min(len(x), max(1, int(math.ceil(2.0 * integrated_autocorrelation_time(x)))))
    seed = int.from_bytes(hashlib.sha256(x.tobytes()).digest()[:8], "little")
    rng = np.random.Generator(np.random.PCG64(seed))
    nblocks = int(math.ceil(len(x) / block))
    starts = rng.integers(0, len(x), size=(replicates, nblocks))
    offsets = np.arange(block)
    indices = (starts[:, :, None] + offsets[None, None, :]) % len(x)
    samples = x[indices.reshape(replicates, -1)[:, :len(x)]]
    means = samples.mean(axis=1)
    alpha = 1.0 - confidence
    return float(np.quantile(means, alpha / 2)), float(np.quantile(means, 1 - alpha / 2)), block


def block_bootstrap_slope_ci(times: Sequence[float], values: Sequence[float], confidence: float = 0.95, replicates: int = 2000) -> tuple[float, float, float, int]:
    """Block-bootstrap OLS slope CI with deterministic analysis-only resampling."""
    t = _finite_array(times, "slope times")
    x = _finite_array(values, "slope values")
    if len(t) != len(x) or len(x) < 4 or np.any(np.diff(t) <= 0):
        raise CriteriaError("slope input must have equal-length, increasing times and >=4 samples")
    tau = integrated_autocorrelation_time(x)
    block = min(len(x), max(1, int(math.ceil(2.0 * tau))))
    seed_data = np.column_stack((t, x)).tobytes()
    seed = int.from_bytes(hashlib.sha256(seed_data).digest()[:8], "little")
    rng = np.random.Generator(np.random.PCG64(seed))
    nblocks = int(math.ceil(len(x) / block))
    starts = rng.integers(0, len(x), size=(replicates, nblocks))
    indices = (starts[:, :, None] + np.arange(block)[None, None, :]) % len(x)
    samples = x[indices.reshape(replicates, -1)[:, :len(x)]]
    slopes = np.empty(replicates)
    centered_t = t - t.mean()
    denom = float(centered_t @ centered_t)
    for idx, row in enumerate(samples):
        slopes[idx] = float(centered_t @ (row - row.mean()) / denom)
    alpha = 1.0 - confidence
    low, high = np.quantile(slopes, [alpha / 2, 1 - alpha / 2])
    point = float(centered_t @ (x - x.mean()) / denom)
    return point, float(low), float(high), block


def multiple_origin_msd(unwrapped_positions_A: Sequence[Sequence[Sequence[float]]], frame_interval_ps: float) -> tuple[np.ndarray, np.ndarray]:
    """Multiple-time-origin MSD; returns lag time (ps) and MSD (A^2)."""
    xyz = _finite_array(unwrapped_positions_A, "unwrapped positions")
    if xyz.ndim != 3 or xyz.shape[2] != 3 or xyz.shape[0] < 4:
        raise CriteriaError("MSD needs >=4 frames of N x 3 unwrapped positions")
    max_lag = xyz.shape[0] // 2
    lag_time = np.arange(1, max_lag + 1, dtype=float) * float(frame_interval_ps)
    msd = np.empty(max_lag, dtype=float)
    for lag in range(1, max_lag + 1):
        delta = xyz[lag:] - xyz[:-lag]
        msd[lag - 1] = float(np.mean(np.sum(delta * delta, axis=2)))
    return lag_time, msd


def rdf_shape_assess(frames: Sequence[FrameMetrics], policy: dict) -> dict:
    """Operationalize the frozen qualitative RDF rule using C-signed limits.

    Required policy keys: minimum_first_peak_fwhm_A,
    secondary_search_start_A, maximum_secondary_peak_prominence_fraction.
    The caller must bind the canonical policy JSON SHA in C's release.
    """
    if not frames:
        raise CriteriaError("RDF shape assessment has no frames")
    required = {
        "minimum_first_peak_fwhm_A",
        "secondary_search_start_A",
        "maximum_secondary_peak_prominence_fraction",
    }
    if not required.issubset(policy):
        raise CriteriaError("C-signed RDF shape policy is incomplete")
    r = np.asarray(frames[0].rdf_r_A, dtype=float)
    curves = np.asarray([f.rdf_g for f in frames], dtype=float)
    if curves.ndim != 2 or curves.shape[1] != len(r) or not np.allclose([f.rdf_r_A for f in frames], r[None, :]):
        raise CriteriaError("RDF grids differ inside a QA window")
    g = curves.mean(axis=0)
    peak_region = np.flatnonzero((r >= 2.6) & (r <= 3.2))
    valley_region = np.flatnonzero((r >= 3.4) & (r <= 4.2))
    if not len(peak_region) or not len(valley_region):
        raise CriteriaError("RDF range does not include the locked first shell")
    peak = int(peak_region[np.argmax(g[peak_region])])
    valley_candidates = valley_region[valley_region > peak]
    if not len(valley_candidates):
        raise CriteriaError("RDF first valley occurs before the first peak")
    valley = int(valley_candidates[np.argmin(g[valley_candidates])])
    half = g[valley] + 0.5 * (g[peak] - g[valley])
    left = peak
    while left > 0 and g[left] >= half:
        left -= 1
    right = peak
    while right + 1 < len(g) and g[right] >= half:
        right += 1
    fwhm = float(r[right] - r[left])
    search_start = float(policy["secondary_search_start_A"])
    secondary = np.flatnonzero((r >= search_start) & (np.arange(len(r)) > valley))
    first_height = float(g[peak] - g[valley])
    secondary_prominences = []
    for idx in secondary:
        if idx <= 0 or idx >= len(g) - 1 or not (g[idx - 1] < g[idx] > g[idx + 1]):
            continue
        left_min = float(np.min(g[max(0, idx - 12):idx + 1]))
        right_min = float(np.min(g[idx:min(len(g), idx + 13)]))
        prominence = max(0.0, float(g[idx] - max(left_min, right_min)))
        secondary_prominences.append(prominence / max(first_height, np.finfo(float).eps))
    max_ratio = max(secondary_prominences, default=0.0)
    fwhm_pass = fwhm >= float(policy["minimum_first_peak_fwhm_A"])
    secondary_pass = max_ratio <= float(policy["maximum_secondary_peak_prominence_fraction"])
    return {
        "first_peak_A": float(r[peak]),
        "first_valley_A": float(r[valley]),
        "first_peak_fwhm_A": fwhm,
        "maximum_secondary_peak_prominence_fraction": float(max_ratio),
        "minimum_first_peak_fwhm_pass": bool(fwhm_pass),
        "secondary_peak_pass": bool(secondary_pass),
        "rdf_shape_pass": bool(fwhm_pass and secondary_pass),
    }


def msd_summary(unwrapped_positions_A: Sequence[Sequence[Sequence[float]]], frame_interval_ps: float) -> dict:
    """4 ps displacement plus a multiple-origin long-lag slope bootstrap."""
    xyz = _finite_array(unwrapped_positions_A, "MSD positions")
    if xyz.ndim != 3 or xyz.shape[2] != 3 or xyz.shape[0] < 8:
        raise CriteriaError("MSD QA requires an adequate unwrapped trajectory interval")
    delta = xyz[-1] - xyz[0]
    delta_msd = float(np.mean(np.sum(delta * delta, axis=1)))
    max_lag = xyz.shape[0] // 2
    n_origins = xyz.shape[0] - max_lag
    per_origin = np.empty((n_origins, max_lag), dtype=float)
    for origin in range(n_origins):
        for lag in range(1, max_lag + 1):
            displacement = xyz[origin + lag] - xyz[origin]
            per_origin[origin, lag - 1] = float(np.mean(np.sum(displacement * displacement, axis=1)))
    lag_steps = np.arange(1, max_lag + 1, dtype=float) * float(frame_interval_ps)
    fit_start = max(1, max_lag // 2) - 1
    fit_t = lag_steps[fit_start:] - lag_steps[fit_start:].mean()
    denom = float(fit_t @ fit_t)
    block = min(n_origins, max(1, int(math.ceil(2.0 * integrated_autocorrelation_time(per_origin[:, -1])))))
    seed = int.from_bytes(hashlib.sha256(per_origin.tobytes()).digest()[:8], "little")
    rng = np.random.Generator(np.random.PCG64(seed))
    reps = 2000
    nblocks = int(math.ceil(n_origins / block))
    starts = rng.integers(0, n_origins, size=(reps, nblocks))
    idx = (starts[:, :, None] + np.arange(block)[None, None, :]) % n_origins
    idx = idx.reshape(reps, -1)[:, :n_origins]
    mean_curve = per_origin.mean(axis=0)
    point = float(fit_t @ (mean_curve[fit_start:] - mean_curve[fit_start:].mean()) / denom)
    slopes = np.empty(reps)
    for k in range(reps):
        curve = per_origin[idx[k]].mean(axis=0)
        slopes[k] = float(fit_t @ (curve[fit_start:] - curve[fit_start:].mean()) / denom)
    low = float(np.quantile(slopes, 0.025))
    high = float(np.quantile(slopes, 0.975))
    return {
        "delta_msd_A2_over_full_interval": delta_msd,
        "multiple_origin_long_lag_slope_A2_per_ps": point,
        "long_lag_slope_ci95_A2_per_ps": [low, high],
        "origin_block_length": block,
        "frame_interval_ps": float(frame_interval_ps),
        "long_lag_fit_range_ps": [float(lag_steps[fit_start]), float(lag_steps[-1])],
        "msd_gate_pass": bool(delta_msd >= 1.0 and low > 0.0),
    }


def stress_drift_gate(frames: Sequence[FrameMetrics], times_ps: Sequence[float]) -> dict:
    if len(frames) != len(times_ps) or len(frames) < 8:
        raise CriteriaError("stress drift gate needs >=8 stress frames with times")
    values = np.asarray([f.stress_eV_A3 for f in frames], dtype=float)
    if values.ndim != 2 or values.shape[1] != 6:
        raise CriteriaError("stress must contain the six ASE Voigt components")
    components = []
    for index in range(6):
        slope, low, high, block = block_bootstrap_slope_ci(times_ps, values[:, index])
        components.append({"component": index, "slope_eV_A3_per_ps": slope, "ci95": [low, high], "block_length": block, "no_resolved_drift": low <= 0 <= high})
    return {"components": components, "stress_drift_pass": all(item["no_resolved_drift"] for item in components)}


def candidate_snapshot_spacing_steps(frames: Sequence[FrameMetrics], thermal_values: Sequence[float], energy_values: Sequence[float]) -> dict:
    """Implement the frozen spacing max(0.5 ps, 5*tau_int) conservatively."""
    if len(frames) < 8:
        raise CriteriaError("candidate spacing requires >=8 post-equilibration analysis frames")
    series = {
        "temperature": np.asarray(thermal_values, dtype=float),
        "potential_energy_per_atom": np.asarray(energy_values, dtype=float),
        "global_q6": np.asarray([f.global_q6 for f in frames], dtype=float),
        "coordination": np.asarray([f.coordination_number for f in frames], dtype=float),
    }
    taus_ps = {}
    for name, values in series.items():
        if len(values) < 4:
            raise CriteriaError(f"not enough {name} values for autocorrelation spacing")
        # All inputs are sampled at 0.05 ps in this plan.
        taus_ps[name] = integrated_autocorrelation_time(values) * 0.05
    tau_max = max(taus_ps.values())
    spacing_ps = max(0.5, 5.0 * tau_max)
    spacing_steps = int(math.ceil(spacing_ps / 0.05) * FRAME_INTERVAL_STEPS)
    return {
        "tau_int_ps_by_observable": taus_ps,
        "tau_int_max_ps": float(tau_max),
        "minimum_snapshot_spacing_ps": float(spacing_ps),
        "snapshot_spacing_steps": spacing_steps,
        "snapshot_spacing_ps_rounded_to_frame_grid": spacing_steps * 0.0005,
        "source_family_policy": "all snapshots remain in this single task family",
    }


def structure_window_summary(frames: Sequence[FrameMetrics], thresholds: dict, rdf_policy: dict | None = None) -> dict:
    if not frames:
        raise CriteriaError("structure window has no frames")
    q6_ok = [f.global_q6 < thresholds["global_Steinhardt_Q6_max"] for f in frames]
    crystal_ok = [
        f.crystalline_fraction <= thresholds["CNA_or_PTM_crystalline_atom_fraction_max_percent"] / 100.0
        and f.largest_crystalline_cluster_atoms <= thresholds["largest_connected_crystalline_cluster_atoms_max"]
        for f in frames
    ]
    q6_fraction = float(np.mean(q6_ok))
    crystal_fraction = float(np.mean(crystal_ok))
    peaks = [f.rdf_first_peak_A for f in frames if f.rdf_first_peak_A is not None]
    valleys = [f.rdf_first_minimum_A for f in frames if f.rdf_first_minimum_A is not None]
    coordination = [f.coordination_number for f in frames if math.isfinite(f.coordination_number)]
    first_peak_range = thresholds["RDF_first_peak_range_A"]
    first_min_range = thresholds["RDF_first_minimum_range_A"]
    coord_range = thresholds["coordination_number_range"]
    rdf_numeric_ok = bool(peaks and valleys and coordination) and (
        all(first_peak_range[0] <= x <= first_peak_range[1] for x in peaks)
        and all(first_min_range[0] <= x <= first_min_range[1] for x in valleys)
        and all(coord_range[0] <= x <= coord_range[1] for x in coordination)
    )
    min_distance = min(f.minimum_pair_distance_A for f in frames)
    shape = rdf_shape_assess(frames, rdf_policy) if rdf_policy is not None else None
    return {
        "frame_count": len(frames),
        "q6_pass_fraction": q6_fraction,
        "crystal_limits_pass_fraction": crystal_fraction,
        "mean_global_q6": float(np.mean([f.global_q6 for f in frames])),
        "mean_crystalline_fraction": float(np.mean([f.crystalline_fraction for f in frames])),
        "maximum_crystalline_fraction": float(np.max([f.crystalline_fraction for f in frames])),
        "maximum_crystalline_cluster_atoms": max(f.largest_crystalline_cluster_atoms for f in frames),
        "rdf_first_peak_mean_A": float(np.mean(peaks)) if peaks else None,
        "rdf_first_minimum_mean_A": float(np.mean(valleys)) if valleys else None,
        "coordination_mean": float(np.mean(coordination)) if coordination else None,
        "minimum_pair_distance_A": float(min_distance),
        "rdf_numeric_ranges_pass": bool(rdf_numeric_ok),
        "rdf_shape_qualitative": "PASS" if shape and shape["rdf_shape_pass"] else ("FAIL" if shape else "REQUIRES_C_OPERATIONAL_POLICY"),
        "rdf_shape": shape,
        "rdf_shape_pass": shape["rdf_shape_pass"] if shape else None,
        "q6_gate_pass": q6_fraction >= thresholds["fraction_frames_meeting_Q6_at_least"],
        "crystal_gate_pass": crystal_fraction >= thresholds["fraction_frames_meeting_crystal_limits_at_least"],
        "minimum_distance_pass": min_distance >= thresholds["minimum_pair_distance_A"],
        "automatic_structure_gate_pass": bool(
            q6_fraction >= thresholds["fraction_frames_meeting_Q6_at_least"]
            and crystal_fraction >= thresholds["fraction_frames_meeting_crystal_limits_at_least"]
            and rdf_numeric_ok
            and min_distance >= thresholds["minimum_pair_distance_A"]
            and shape is not None
            and shape["rdf_shape_pass"]
        ) if shape else None,
        "note": "C frozen lock specifies a broad first peak and no persistent sharp fcc/hcp multishell order but gives no numerical decision rule; automatic pass remains blocked until C signs an operational rule.",
    }


def paired_window_thermal_gate(first: dict, second: dict, thresholds: dict, target_K: float) -> dict:
    """Apply C's numerical thermal/energy gates to adjacent 2 ps windows."""
    t1 = _finite_array(first["temperature_K"], "first-window temperatures")
    t2 = _finite_array(second["temperature_K"], "second-window temperatures")
    e1 = _finite_array(first["potential_energy_eV_per_atom"], "first-window energies")
    e2 = _finite_array(second["potential_energy_eV_per_atom"], "second-window energies")
    times = np.concatenate((
        np.asarray(first["time_ps"], dtype=float),
        np.asarray(second["time_ps"], dtype=float),
    ))
    temp = np.concatenate((t1, t2))
    mean1, mean2 = float(t1.mean()), float(t2.mean())
    tolerance = thresholds["temperature_mean_tolerance_percent"] / 100.0
    adjacent_tol = thresholds["adjacent_window_mean_difference_tolerance_percent"] / 100.0
    slope, slope_low, slope_high, temp_block = block_bootstrap_slope_ci(times, temp)
    ci1 = block_bootstrap_mean_ci(t1)
    ci2 = block_bootstrap_mean_ci(t2)
    energy_delta = abs(float(e2.mean() - e1.mean()))
    energy_limit = thresholds["potential_energy_per_atom_adjacent_window_difference_limit_eV"]
    return {
        "temperature_mean_first_K": mean1,
        "temperature_mean_second_K": mean2,
        "temperature_std_first_K": float(t1.std(ddof=1)),
        "temperature_std_second_K": float(t2.std(ddof=1)),
        "temperature_ci95_first_K": list(ci1[:2]),
        "temperature_ci95_second_K": list(ci2[:2]),
        "temperature_block_length_samples": max(ci1[2], ci2[2], temp_block),
        "temperature_slope_K_per_ps": slope,
        "temperature_slope_ci95_K_per_ps": [slope_low, slope_high],
        "potential_energy_mean_first_eV_per_atom": float(e1.mean()),
        "potential_energy_mean_second_eV_per_atom": float(e2.mean()),
        "potential_energy_adjacent_delta_eV_per_atom": energy_delta,
        "temperature_means_in_range": abs(mean1 - target_K) <= tolerance * target_K and abs(mean2 - target_K) <= tolerance * target_K,
        "adjacent_temperature_difference_in_range": abs(mean2 - mean1) <= adjacent_tol * target_K,
        "temperature_slope_in_range": abs(slope) <= thresholds["absolute_temperature_trend_limit_K_per_ps"],
        "temperature_ci_covers_target": ci1[0] <= target_K <= ci1[1] and ci2[0] <= target_K <= ci2[1],
        "energy_adjacent_difference_in_range": energy_delta <= energy_limit,
        "temperature_gate_pass": (
            abs(mean1 - target_K) <= tolerance * target_K
            and abs(mean2 - target_K) <= tolerance * target_K
            and abs(mean2 - mean1) <= adjacent_tol * target_K
            and abs(slope) <= thresholds["absolute_temperature_trend_limit_K_per_ps"]
            and ci1[0] <= target_K <= ci1[1]
            and ci2[0] <= target_K <= ci2[1]
        ),
        "energy_gate_pass": energy_delta <= energy_limit,
    }


def monitor_record_to_json(metrics: FrameMetrics) -> dict:
    return asdict(metrics)
