#!/usr/bin/env python3
"""C QA V2 liquid checks for the locked finite Ag64/rho100 cell.

PTM is delegated to OVITO's PolyhedralTemplateMatchingModifier.  This module
never substitutes CNA when PTM is absent and all unresolved inputs fail closed.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Sequence

import numpy as np
from ase import Atoms
from ase.io import write
from ase.neighborlist import neighbor_list
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths
from scipy.special import sph_harm_y


PACKAGE_DIR = Path(__file__).resolve().parent
SCIENCE_LOCK_PATH = PACKAGE_DIR / "LIQUID_QA_SCIENCE_LOCK_V2_1.json"
SCIENCE_LOCK_SHA256 = "848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2"
LOCK = json.loads(SCIENCE_LOCK_PATH.read_text(encoding="utf-8"))
FRAME_INTERVAL_PS = 0.05
STEP_INTERVAL = 100
BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_SEED = 20261012
ASE_STRESS_TO_GPA = 160.21766208
STRESS_SAMPLE_INTERVAL_PS = 0.01
STRESS_MIN_STATIONARY_DURATION_PS = 4.0


class CriteriaError(RuntimeError):
    pass


@dataclass(frozen=True)
class FrameMetrics:
    step: int
    time_ps: float
    temperature_K: float
    potential_energy_eV_per_atom: float
    global_q6: float
    mean_local_q6: float
    ptm_status: str
    ptm_crystalline_fraction: float | None
    ptm_largest_cluster_atoms: int | None
    ptm_fcc_fraction: float | None
    ptm_hcp_fraction: float | None
    ptm_mean_rmsd: float | None
    ptm_version: str | None
    ptm_rmsd_cutoff: float
    pbc_verified: bool
    coordination_number: float | None
    minimum_pair_distance_A: float
    rdf_counts: tuple[int, ...]
    stress_eV_A3: tuple[float, ...]


def _finite_array(value, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if not np.isfinite(arr).all():
        raise CriteriaError(f"{name} contains non-finite values")
    return arr


def _locked_edges(lock: dict = LOCK) -> np.ndarray:
    rdf = lock["rdf"]
    text = rdf["histogram"]["edges_angstrom"]
    rmax = float(lock["finite_cell"]["rdf_rmax_angstrom"])
    if text != "0.00 through 5.35 inclusive, fixed 0.05 A bins":
        raise CriteriaError("unrecognized V2 RDF edge specification")
    return np.linspace(0.0, rmax, int(round(rmax / 0.05)) + 1)


def validate_locked_atoms(atoms: Atoms, lock: dict = LOCK) -> None:
    n = int(lock["finite_cell"]["atoms"])
    locked_lengths = np.asarray(lock["finite_cell"]["cell_angstrom"], dtype=float)
    cell = _finite_array(np.asarray(atoms.cell), "cell")
    pos = _finite_array(atoms.get_positions(), "positions")
    if len(atoms) != n or atoms.get_chemical_formula() != "Ag64" or pos.shape != (n, 3):
        raise CriteriaError("locked composition/position shape mismatch")
    if not np.array_equal(np.asarray(atoms.pbc, dtype=bool), np.asarray(lock["finite_cell"]["pbc"], dtype=bool)):
        raise CriteriaError("locked PBC identity mismatch")
    if locked_lengths.shape != (3,):
        raise CriteriaError("locked cell lengths must contain three values")
    expected = np.diag(locked_lengths)
    if cell.shape != (3, 3) or not np.allclose(cell, expected, rtol=0.0, atol=1e-12):
        raise CriteriaError("locked cell identity mismatch")
    rmax = float(lock["finite_cell"]["rdf_rmax_angstrom"])
    shortest = float(np.min(np.linalg.norm(cell, axis=1)))
    if not rmax < shortest / 2.0:
        raise CriteriaError("RDF radial maximum is not strictly inside the complete-shell radius")


def rdf_counts(atoms: Atoms, lock: dict = LOCK) -> tuple[np.ndarray, float]:
    """Ordered-pair histogram with minimum-image distances under full PBC."""
    validate_locked_atoms(atoms, lock)
    distances = _finite_array(atoms.get_all_distances(mic=True), "minimum-image distances")
    n = len(atoms)
    ordered = distances[~np.eye(n, dtype=bool)]
    edges = _locked_edges(lock)
    rmax = float(edges[-1])
    counts = np.histogram(ordered[ordered <= rmax], bins=edges)[0].astype(np.int64)
    minimum = float(np.min(ordered))
    return counts, minimum


def normalized_rdf(counts: Sequence[int], frame_count: int, lock: dict = LOCK) -> tuple[np.ndarray, np.ndarray]:
    edges = _locked_edges(lock)
    hist = np.asarray(counts, dtype=float)
    if frame_count < 1 or hist.shape != (len(edges) - 1,) or not np.isfinite(hist).all() or np.any(hist < 0):
        raise CriteriaError("invalid RDF counts or frame count")
    n = int(lock["finite_cell"]["atoms"])
    volume = float(lock["finite_cell"]["volume_angstrom3"])
    shells = (4.0 * math.pi / 3.0) * (edges[1:] ** 3 - edges[:-1] ** 3)
    denom = frame_count * n * (n - 1) * shells / volume
    g = hist / denom
    if not np.isfinite(g).all():
        raise CriteriaError("normalized RDF is non-finite")
    return (edges[:-1] + edges[1:]) / 2.0, g


def local_q6(atoms: Atoms, cutoff_A: float = 3.8) -> tuple[float, np.ndarray]:
    """V2 mean local q6: normalize each atom's neighbor-averaged Y6m first."""
    validate_locked_atoms(atoms)
    i, _j, shifts, distances = neighbor_list("ijSd", atoms, cutoff_A, self_interaction=False)
    if len(i) == 0 or not np.isfinite(distances).all() or np.any(distances <= 0):
        raise CriteriaError("Q6 has no valid periodic neighbors")
    cell = np.asarray(atoms.cell, dtype=float)
    vectors = np.asarray(atoms.positions)[_j] + np.asarray(shifts) @ cell - np.asarray(atoms.positions)[i]
    theta = np.arccos(np.clip(vectors[:, 2] / distances, -1.0, 1.0))
    phi = np.arctan2(vectors[:, 1], vectors[:, 0])
    degrees = np.bincount(i, minlength=len(atoms))
    if np.any(degrees == 0):
        raise CriteriaError("Q6 found an atom with no neighbors inside 3.8 A")
    qlm = np.empty((len(atoms), 13), dtype=np.complex128)
    for column, m in enumerate(range(-6, 7)):
        y = sph_harm_y(6, m, theta, phi)
        qlm[:, column] = np.bincount(i, weights=y.real, minlength=len(atoms)) + 1j * np.bincount(i, weights=y.imag, minlength=len(atoms))
    qlm /= degrees[:, None]
    q6_i = np.sqrt((4.0 * math.pi / 13.0) * np.sum(np.abs(qlm) ** 2, axis=1))
    _finite_array(q6_i, "local q6")
    return float(q6_i.mean()), q6_i


def global_q6(atoms: Atoms, cutoff_A: float = 3.8) -> float:
    """Legacy coherent Steinhardt Q6; must retain the archived reduction order."""
    validate_locked_atoms(atoms)
    i, j, shifts, distances = neighbor_list("ijSd", atoms, cutoff_A, self_interaction=False)
    if len(i) == 0 or not np.isfinite(distances).all() or np.any(distances <= 0):
        raise CriteriaError("global Q6 has no valid periodic neighbors")
    vectors = np.asarray(atoms.positions)[j] + np.asarray(shifts) @ np.asarray(atoms.cell) - np.asarray(atoms.positions)[i]
    theta = np.arccos(np.clip(vectors[:, 2] / distances, -1.0, 1.0))
    phi = np.arctan2(vectors[:, 1], vectors[:, 0])
    degree = np.bincount(i, minlength=len(atoms))
    if np.any(degree == 0):
        raise CriteriaError("global Q6 found an atom with no neighbors inside 3.8 A")
    qlm = np.zeros((len(atoms), 13), dtype=np.complex128)
    for column, m in enumerate(range(-6, 7)):
        y = sph_harm_y(6, m, theta, phi)
        qlm[:, column] = np.bincount(i, weights=y.real, minlength=len(atoms)) + 1j * np.bincount(i, weights=y.imag, minlength=len(atoms))
    qlm /= degree[:, None]
    mean_qlm = qlm.mean(axis=0)
    value = math.sqrt(4.0 * math.pi / 13.0 * float(np.sum(np.abs(mean_qlm) ** 2)))
    if not math.isfinite(value):
        raise CriteriaError("global Q6 is non-finite")
    return value


def _largest_periodic_cluster(atoms: Atoms, crystalline: np.ndarray, cutoff_A: float = 3.8) -> int:
    n = len(atoms)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    ii, jj = neighbor_list("ij", atoms, cutoff_A, self_interaction=False)
    for a, b in zip(ii, jj):
        if crystalline[int(a)] and crystalline[int(b)]:
            union(int(a), int(b))
    sizes: dict[int, int] = {}
    for i, flag in enumerate(crystalline):
        if flag:
            root = find(i)
            sizes[root] = sizes.get(root, 0) + 1
    return max(sizes.values(), default=0)


def ptm_classify(atoms: Atoms, lock: dict = LOCK) -> dict:
    """Run genuine OVITO PTM with the C-locked RMSD cutoff and periodic cell."""
    classifier = lock["order"]["other_crystalline_order_gates"]["classifier"]
    match = re.search(r"RMSD cutoff ([0-9]+(?:\.[0-9]+)?)", classifier)
    if not match:
        return {"status": "PTM_UNAVAILABLE", "reason": "locked PTM RMSD cutoff is missing"}
    cutoff = float(match.group(1))
    try:
        if not np.array_equal(np.asarray(atoms.pbc, dtype=bool), np.asarray([True, True, True])):
            return {"status": "PERIODICITY_FAIL", "reason": "PTM input is not fully periodic"}
        import importlib.metadata
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from ovito.io import import_file
        from ovito.modifiers import PolyhedralTemplateMatchingModifier
        ovito_version = importlib.metadata.version("ovito")
    except Exception as exc:
        return {"status": "PTM_UNAVAILABLE", "reason": f"OVITO PTM import failed: {type(exc).__name__}: {exc}"}
    try:
        with tempfile.TemporaryDirectory(prefix="ag64-ptm-frame-") as temp:
            xyz = Path(temp) / "frame.extxyz"
            write(xyz, atoms, format="extxyz")
            pipeline = import_file(str(xyz))
            modifier = PolyhedralTemplateMatchingModifier()
            modifier.rmsd_cutoff = cutoff
            modifier.output_rmsd = True
            pipeline.modifiers.append(modifier)
            data = pipeline.compute()
            if tuple(data.cell.pbc) != tuple(bool(x) for x in atoms.pbc):
                return {"status": "PERIODICITY_FAIL", "reason": "OVITO did not preserve the input PBC"}
            # ASE stores cell vectors as rows; OVITO stores them as columns.
            ovito_cell = np.asarray(data.cell.matrix, dtype=float)[:, :3].T
            if not np.allclose(ovito_cell, np.asarray(atoms.cell), rtol=0.0, atol=1e-9):
                return {"status": "PERIODICITY_FAIL", "reason": "OVITO did not preserve the input cell"}
            props = data.particles
            if "Structure Type" not in props or "RMSD" not in props:
                return {"status": "PTM_UNAVAILABLE", "reason": "OVITO PTM omitted required Structure Type or RMSD output"}
            types = np.asarray(props["Structure Type"].array, dtype=int)
            rmsd = _finite_array(props["RMSD"].array, "PTM RMSD")
            if types.shape != (len(atoms),) or rmsd.shape != (len(atoms),):
                return {"status": "PTM_UNAVAILABLE", "reason": "OVITO PTM output atom count differs"}
            fcc_id = int(PolyhedralTemplateMatchingModifier.Type.FCC.value)
            hcp_id = int(PolyhedralTemplateMatchingModifier.Type.HCP.value)
            fcc = types == fcc_id
            hcp = types == hcp_id
            crystalline = (fcc | hcp) & (rmsd <= cutoff + 1e-12)
            return {
                "status": "PASS",
                "version": ovito_version,
                "implementation": "OVITO PolyhedralTemplateMatchingModifier",
                "rmsd_cutoff": cutoff,
                "fcc_fraction": float(np.count_nonzero(fcc & crystalline) / len(atoms)),
                "hcp_fraction": float(np.count_nonzero(hcp & crystalline) / len(atoms)),
                "crystalline_fraction": float(np.count_nonzero(crystalline) / len(atoms)),
                "largest_cluster_atoms": _largest_periodic_cluster(atoms, crystalline, 3.8),
                "mean_rmsd": float(rmsd.mean()),
                "pbc_verified": True,
            }
    except Exception as exc:
        return {"status": "PTM_ERROR", "reason": f"OVITO PTM compute failed: {type(exc).__name__}: {exc}"}


def frame_metrics(atoms: Atoms, step: int, temperature_K: float, potential_energy_eV: float,
                  stress_eV_A3: Sequence[float] | None = None, lock: dict = LOCK) -> FrameMetrics:
    validate_locked_atoms(atoms, lock)
    temp = float(temperature_K)
    energy_per_atom = float(potential_energy_eV) / len(atoms)
    if not np.isfinite([temp, energy_per_atom]).all():
        raise CriteriaError("frame temperature or energy is non-finite")
    q6_global = global_q6(atoms)
    q6_mean, _q6_i = local_q6(atoms)
    counts, mind = rdf_counts(atoms, lock)
    ptm = ptm_classify(atoms, lock)
    stress = tuple(float(x) for x in (stress_eV_A3 if stress_eV_A3 is not None else [0.0] * 6))
    _finite_array(stress, "stress")
    cn = None
    # A frame-level CN is diagnostic; the pooled phase valley is used for the
    # actual gate and is recomputed inside analyze_rdf_phase().
    point = _rdf_point(counts[None, :], lock)
    if point["status"] == "PASS":
        cn = float(point["coordination_number"])
    return FrameMetrics(
        step=int(step), time_ps=int(step) * 0.0005, temperature_K=temp,
        potential_energy_eV_per_atom=energy_per_atom, global_q6=q6_global, mean_local_q6=q6_mean,
        ptm_status=ptm["status"], ptm_crystalline_fraction=ptm.get("crystalline_fraction"),
        ptm_largest_cluster_atoms=ptm.get("largest_cluster_atoms"),
        ptm_fcc_fraction=ptm.get("fcc_fraction"), ptm_hcp_fraction=ptm.get("hcp_fraction"),
        ptm_mean_rmsd=ptm.get("mean_rmsd"), ptm_version=ptm.get("version"),
        ptm_rmsd_cutoff=0.10, pbc_verified=bool(ptm.get("pbc_verified", False)),
        coordination_number=cn, minimum_pair_distance_A=mind,
        rdf_counts=tuple(int(x) for x in counts), stress_eV_A3=stress,
    )


def monitor_record_to_json(metrics: FrameMetrics) -> dict:
    return asdict(metrics)


def frame_metrics_from_json(row: dict) -> FrameMetrics:
    data = dict(row)
    data["rdf_counts"] = tuple(int(x) for x in data["rdf_counts"])
    data["stress_eV_A3"] = tuple(float(x) for x in data["stress_eV_A3"])
    return FrameMetrics(**data)


def _rdf_smooth(g: np.ndarray, lock: dict = LOCK) -> np.ndarray:
    if lock["rdf"]["histogram"]["smoothing_for_extrema_only"] != "Gaussian filter sigma=0.10 A (2 bins); raw unsmoothed g(r) is used for coordination integration":
        raise CriteriaError("unrecognized RDF smoothing lock")
    return gaussian_filter1d(g, sigma=2.0, mode="nearest")


def _peak_valley(r: np.ndarray, raw_g: np.ndarray, lock: dict = LOCK) -> dict:
    smooth = _rdf_smooth(raw_g, lock)
    p_lo, p_hi = map(float, lock["rdf"]["peak_valley"]["peak_search_range_angstrom"])
    v_lo, v_hi = map(float, lock["rdf"]["peak_valley"]["valley_search_range_angstrom"])
    peak_candidates, _ = find_peaks(smooth)
    peak_candidates = [int(k) for k in peak_candidates if p_lo <= r[k] <= p_hi]
    if not peak_candidates:
        return {"status": "RDF_UNRESOLVED", "reason": "no local maximum in the locked first-peak range", "smoothed_g": smooth}
    peak_idx = max(peak_candidates, key=lambda k: smooth[k])
    valley_candidates, _ = find_peaks(-smooth)
    valley_candidates = [int(k) for k in valley_candidates if r[k] >= max(v_lo, r[peak_idx]) and r[k] <= v_hi]
    if not valley_candidates:
        return {"status": "RDF_UNRESOLVED", "reason": "no local minimum after first peak in the locked valley range", "smoothed_g": smooth}
    valley_idx = min(valley_candidates, key=lambda k: smooth[k])
    return {
        "status": "PASS", "peak_index": peak_idx, "valley_index": valley_idx,
        "peak_A": float(r[peak_idx]), "valley_A": float(r[valley_idx]),
        "peak_height_raw": float(raw_g[peak_idx]), "valley_height_raw": float(raw_g[valley_idx]),
        "peak_valley_delta_raw": float(raw_g[peak_idx] - raw_g[valley_idx]),
        "smoothed_g": smooth,
    }


def _coordination_number(r: np.ndarray, raw_g: np.ndarray, valley_A: float, lock: dict = LOCK) -> float:
    n = int(lock["finite_cell"]["atoms"])
    volume = float(lock["finite_cell"]["volume_angstrom3"])
    inside = r < valley_A
    rr = np.concatenate(([0.0], r[inside], [valley_A]))
    gg = np.concatenate(([0.0], raw_g[inside], [float(np.interp(valley_A, r, raw_g))]))
    integral = float(np.trapezoid(gg * rr**2, rr))
    return float(4.0 * math.pi * (n - 1) / volume * integral)


def _rdf_point(count_rows: np.ndarray, lock: dict = LOCK) -> dict:
    rows = np.asarray(count_rows, dtype=float)
    if rows.ndim != 2 or len(rows) == 0:
        return {"status": "RDF_UNRESOLVED", "reason": "no frame histograms"}
    r, g = normalized_rdf(rows.sum(axis=0), len(rows), lock)
    pv = _peak_valley(r, g, lock)
    if pv["status"] != "PASS":
        return {k: v for k, v in pv.items() if k != "smoothed_g"}
    pv.pop("smoothed_g", None)
    pv["coordination_number"] = _coordination_number(r, g, pv["valley_A"], lock)
    pv["r_A"], pv["raw_g"], pv["smoothed_g"] = r, g, gaussian_filter1d(g, 2.0, mode="nearest")
    return pv


def integrated_autocorrelation_time(values: Sequence[float], dt_ps: float = FRAME_INTERVAL_PS) -> float:
    x = _finite_array(values, "autocorrelation series")
    if x.ndim != 1 or len(x) < 3:
        raise CriteriaError("autocorrelation needs at least three finite samples")
    centered = x - float(np.mean(x))
    variance = float(np.dot(centered, centered) / len(centered))
    if variance <= np.finfo(float).eps:
        return 0.5 * dt_ps
    total = 0.0
    for lag in range(1, len(x)):
        rho = float(np.dot(centered[:-lag], centered[lag:]) / (len(x) - lag) / variance)
        if rho <= 0.0:
            break
        total += rho
    return float(dt_ps * (0.5 + total))


def _circular_block_indices(n: int, block_frames: int, replicates: int,
                            rng: np.random.Generator) -> np.ndarray:
    if n < 1 or block_frames < 1 or block_frames > n:
        raise CriteriaError("invalid circular moving-block dimensions")
    blocks_needed = int(math.ceil(n / block_frames))
    starts = rng.integers(0, n, size=(replicates, blocks_needed))
    offsets = np.arange(block_frames)
    indices = (starts[:, :, None] + offsets[None, None, :]) % n
    return indices.reshape(replicates, -1)[:, :n]


def _stress_series_seed(times_ps: np.ndarray, values: np.ndarray, name: str) -> int:
    payload = name.encode("utf-8") + b"\0" + np.column_stack((times_ps, values)).astype("<f8").tobytes()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _stress_observable_summary(times_ps: np.ndarray, values: np.ndarray, name: str,
                               dt_ps: float, replicates: int) -> dict:
    n = len(values)
    centered_t = times_ps - float(np.mean(times_ps))
    denominator = float(np.dot(centered_t, centered_t))
    if denominator <= 0.0:
        return {"status": "QA_HOLD", "reason": "stress trend time axis is degenerate", "sample_count": n}
    slope = float(np.dot(centered_t, values - float(np.mean(values))) / denominator)
    fitted = float(np.mean(values)) + slope * centered_t
    residual = values - fitted
    residual -= float(np.mean(residual))
    # Estimate correlation in the detrended residual series so a real trend is
    # not mistaken for a near-unit-root autocorrelation and made untestable.
    tau_ps = integrated_autocorrelation_time(residual, dt_ps)
    block_frames = max(1, int(math.ceil(2.0 * tau_ps / dt_ps)))
    effective_blocks = n // block_frames
    if effective_blocks < 8:
        return {
            "status": "PENDING", "reason": "fewer than eight complete autocorrelation-sized blocks",
            "sample_count": n, "mean": float(np.mean(values)),
            "mean_ci95": None, "mean_block_uncertainty_95pct_half_width": None,
            "slope_per_ps": slope, "slope_ci95_per_ps": None,
            "statistically_resolved_drift": None,
            "tau_int_ps": tau_ps, "block_length_samples": block_frames,
            "complete_effective_blocks": effective_blocks,
        }

    rng = np.random.Generator(np.random.PCG64(_stress_series_seed(times_ps, values, name)))
    indices = _circular_block_indices(n, block_frames, replicates, rng)
    mean_replicates = values[indices].mean(axis=1)
    mean_ci = np.quantile(mean_replicates, [0.025, 0.975])
    slope_replicates = np.empty(replicates, dtype=float)
    if np.all(values == values[0]):
        # Avoid classifying floating-point cancellation noise as a trend in an
        # exactly constant synthetic/recorded series.
        slope = 0.0
        slope_ci = np.asarray([0.0, 0.0])
    else:
        for replicate, sample_indices in enumerate(indices):
            bootstrap_values = fitted + residual[sample_indices]
            slope_replicates[replicate] = float(
                np.dot(centered_t, bootstrap_values - float(np.mean(bootstrap_values))) / denominator
            )
        slope_ci = np.quantile(slope_replicates, [0.025, 0.975])
    no_resolved_drift = bool(slope_ci[0] <= 0.0 <= slope_ci[1])
    mean = float(np.mean(values))
    return {
        "status": "PASS" if no_resolved_drift else "FAIL_DRIFT",
        "sample_count": n,
        "mean": mean,
        "mean_ci95": [float(mean_ci[0]), float(mean_ci[1])],
        "mean_block_uncertainty_95pct_half_width": float((mean_ci[1] - mean_ci[0]) / 2.0),
        "slope_per_ps": slope,
        "slope_ci95_per_ps": [float(slope_ci[0]), float(slope_ci[1])],
        "statistically_resolved_drift": not no_resolved_drift,
        "significance_rule": "two-sided 95% residual circular moving-block-bootstrap slope interval excludes zero",
        "tau_int_ps": tau_ps,
        "block_length_samples": block_frames,
        "complete_effective_blocks": effective_blocks,
        "bootstrap_replicates": replicates,
        "bootstrap_seed": _stress_series_seed(times_ps, values, name),
    }


def pressure_stress_gate(times_ps: Sequence[float] | None,
                         stress_eV_A3: Sequence[Sequence[float]] | None,
                         lock: dict = LOCK) -> dict:
    """Require stationary pressure and all six stress components to have no resolved drift.

    ASE stress is tensile-positive in eV/Angstrom^3; hydrostatic pressure is
    therefore P = -trace(stress)/3 (positive for compression).  The gate has
    no absolute pressure acceptance threshold.
    """
    if times_ps is None or stress_eV_A3 is None:
        return {"status": "QA_HOLD", "pressure_stress_gate_pass": False,
                "reason": "pressure/stress input missing"}
    try:
        times = _finite_array(times_ps, "stress times")
        stress = _finite_array(stress_eV_A3, "six-component stress")
    except (CriteriaError, TypeError, ValueError) as exc:
        return {"status": "QA_HOLD", "pressure_stress_gate_pass": False,
                "reason": f"invalid pressure/stress input: {exc}"}
    if times.ndim != 1 or stress.ndim != 2 or stress.shape != (len(times), 6):
        return {"status": "QA_HOLD", "pressure_stress_gate_pass": False,
                "reason": "stress input must contain six finite Voigt components per time"}
    n = len(times)
    if n == 0 and stress.size == 0:
        return {"status": "PENDING", "pressure_stress_gate_pass": False,
                "reason": "no eligible fixed-temperature stress samples", "sample_count": 0}
    if n < 2 or np.any(np.diff(times) <= 0.0):
        return {"status": "PENDING", "pressure_stress_gate_pass": False,
                "reason": "not enough increasing stress samples", "sample_count": n}
    differences = np.diff(times)
    dt_ps = float(np.median(differences))
    if not math.isfinite(dt_ps) or dt_ps <= 0.0 or not np.allclose(
        differences, dt_ps, rtol=1e-5, atol=1e-9
    ):
        return {"status": "QA_HOLD", "pressure_stress_gate_pass": False,
                "reason": "stress samples are not on a regular time grid", "sample_count": n}
    if not math.isclose(dt_ps, STRESS_SAMPLE_INTERVAL_PS, rel_tol=0.0, abs_tol=1e-9):
        return {"status": "QA_HOLD", "pressure_stress_gate_pass": False,
                "reason": "stress sample cadence differs from the frozen 0.01 ps output interval",
                "sample_count": n, "sample_interval_ps": dt_ps}
    duration_ps = float(times[-1] - times[0] + dt_ps)
    if duration_ps + 1e-12 < STRESS_MIN_STATIONARY_DURATION_PS:
        return {"status": "PENDING", "pressure_stress_gate_pass": False,
                "reason": "less than the complete 4 ps fixed-temperature QA interval",
                "sample_count": n, "duration_ps": duration_ps}

    diagonal = stress[:, :3]
    pressure = -np.sum(diagonal, axis=1) / 3.0
    names = ("pressure", "stress_xx", "stress_yy", "stress_zz", "stress_yz", "stress_xz", "stress_xy")
    series = (pressure,) + tuple(stress[:, index] for index in range(6))
    reps = int(lock["rdf"]["bootstrap"]["replicates"])
    summaries = {
        name: _stress_observable_summary(times, values, name, dt_ps, reps)
        for name, values in zip(names, series)
    }
    pending = [name for name, item in summaries.items() if item["status"] == "PENDING"]
    invalid = [name for name, item in summaries.items() if item["status"] == "QA_HOLD"]
    drift = [name for name, item in summaries.items() if item["status"] == "FAIL_DRIFT"]
    if invalid:
        status, reason = "QA_HOLD", "one or more stress statistics are not reliable"
    elif drift:
        status, reason = "FAIL_DRIFT", "statistically resolved pressure/stress drift"
    elif pending:
        status, reason = "PENDING", "one or more stress series have too few effective blocks"
    else:
        status, reason = "PASS", None
    return {
        "status": status,
        "pressure_stress_gate_pass": status == "PASS",
        "reason": reason,
        "pressure_definition": "P=-(sigma_xx+sigma_yy+sigma_zz)/3",
        "stress_convention": "ASE positive tensile; pressure positive compressive",
        "stress_unit": "eV/Angstrom^3",
        "pressure_unit": "eV/Angstrom^3",
        "pressure_mean_GPa": float(np.mean(pressure) * ASE_STRESS_TO_GPA),
        "duration_ps": duration_ps,
        "sample_interval_ps": dt_ps,
        "sample_count": n,
        "bootstrap_method": "circular moving-block bootstrap; OLS detrended residual blocks for slope CI",
        "bootstrap_replicates": reps,
        "confidence_level": 0.95,
        "multiple_comparison_policy": "each of pressure and six components uses its two-sided 95% interval; any resolved drift fails closed",
        "confirmation_required": "C review of slope bootstrap and significance rule before release",
        "series": summaries,
    }


def _slope(values: np.ndarray, dt_ps: float) -> float:
    y = np.asarray(values, dtype=float)
    x = np.arange(len(y), dtype=float) * dt_ps
    xc = x - x.mean()
    denom = float(np.dot(xc, xc))
    return float(np.dot(xc, y - y.mean()) / denom) if denom > 0 else math.nan


def _moving_block_indices(n: int, block_frames: int, replicates: int, rng: np.random.Generator) -> np.ndarray:
    if block_frames < 1 or n < block_frames:
        raise CriteriaError("block length exceeds the available data")
    blocks_needed = int(math.ceil(n / block_frames))
    result = np.empty((replicates, n), dtype=np.int32)
    max_start = n - block_frames
    for b in range(replicates):
        starts = rng.integers(0, max_start + 1, size=blocks_needed)
        sequence = np.concatenate([np.arange(s, s + block_frames) for s in starts])
        result[b] = sequence[:n]
    return result


def _bootstrap_ci(series: np.ndarray, indices: np.ndarray, kind: str = "mean", dt_ps: float = FRAME_INTERVAL_PS) -> list[float]:
    series = _finite_array(series, "bootstrap series")
    vals = np.empty(len(indices), dtype=float)
    if kind == "slope":
        # Resample contiguous residual blocks while keeping the original time
        # axis fixed; regressing a block-reordered series would create spurious
        # trends from the arbitrary ordering of resampled blocks.
        x = np.arange(len(series), dtype=float) * dt_ps
        slope = _slope(series, dt_ps)
        fitted = float(np.mean(series)) + slope * (x - float(np.mean(x)))
        residual = series - fitted
        residual -= float(np.mean(residual))
        for i, idx in enumerate(indices):
            vals[i] = _slope(fitted + residual[idx], dt_ps)
    else:
        for i, idx in enumerate(indices):
            vals[i] = float(np.mean(series[idx]))
    return [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))]


def _window_point(frames: Sequence[FrameMetrics], bounds: Sequence[float], lock: dict = LOCK) -> dict:
    start, end = map(float, bounds)
    selected = [f for f in frames if start <= f.time_ps < end]
    if len(selected) < int(round((end - start) / FRAME_INTERVAL_PS)) - 1:
        return {"status": "QA_HOLD", "reason": "final 2 ps window has missing frames", "frame_count": len(selected)}
    counts = np.asarray([f.rdf_counts for f in selected], dtype=np.int64)
    point = _rdf_point(counts, lock)
    if point.get("status") != "PASS":
        return {"status": "RDF_UNRESOLVED", "reason": point.get("reason"), "frame_count": len(selected)}
    peak_lo, peak_hi = lock["rdf"]["peak_valley"]["accepted_peak_range_angstrom"]
    valley_lo, valley_hi = lock["rdf"]["peak_valley"]["accepted_valley_range_angstrom"]
    cn_lo, cn_hi = lock["rdf"]["coordination_number"]["accepted_range"]
    point["frame_count"] = len(selected)
    point["point_gate_pass"] = bool(
        peak_lo <= point["peak_A"] <= peak_hi
        and valley_lo <= point["valley_A"] <= valley_hi
        and cn_lo <= point["coordination_number"] <= cn_hi
    )
    point["status"] = "PASS" if point["point_gate_pass"] else "FAIL"
    return point


def q6_regression_report(global_q6_frames: Sequence[float], local_q6_frames: Sequence[float] | None = None,
                         lock: dict = LOCK) -> dict:
    values = _finite_array(global_q6_frames, "legacy global-Q6 regression")
    regression = lock["order"]["legacy_global_q6"]["historical_regression"]
    target = float(regression["expected_mean"])
    tolerance = float(regression["tolerance"])
    mean = float(values.mean()) if len(values) else math.nan
    local_mean = float(np.mean(_finite_array(local_q6_frames, "local Q6 diagnostic"))) if local_q6_frames is not None and len(local_q6_frames) else None
    expected_count = int(regression["frame_count"])
    passed = bool(len(values) == expected_count and abs(mean - target) <= tolerance)
    return {
        "frame_count": int(len(values)), "legacy_global_q6_mean": mean,
        "local_q6_mean_report_only": local_mean,
        "documented_target": target, "tolerance": tolerance,
        "regression_pass": passed,
        "status": "PASS" if passed else "LEGACY_GLOBAL_Q6_REGRESSION_HOLD",
        "interpretation": "Only legacy global Q6 is compared to the historical global-Q6 value; local Q6 is diagnostic only.",
    }


def analyze_rdf_phase(frames: Sequence[FrameMetrics], final_windows_ps: Sequence[Sequence[float]],
                      lock: dict = LOCK) -> dict:
    if not frames:
        return {"status": "QA_HOLD", "reason": "no stationary phase frames"}
    if any(f.ptm_status != "PASS" or not f.pbc_verified for f in frames):
        bad = sorted({f.ptm_status for f in frames if f.ptm_status != "PASS" or not f.pbc_verified})
        return {"status": "QA_HOLD", "reason": "PTM unavailable/error or periodicity unverified", "ptm_statuses": bad}
    fields = ["temperature_K", "potential_energy_eV_per_atom", "mean_local_q6", "ptm_crystalline_fraction"]
    try:
        series = {name: _finite_array([getattr(f, name) for f in frames], name) for name in fields}
        cn_point = _rdf_point(np.asarray([f.rdf_counts for f in frames]), lock)
    except (CriteriaError, TypeError, ValueError) as exc:
        return {"status": "QA_HOLD", "reason": f"required metric missing/non-finite: {exc}"}
    if cn_point.get("status") != "PASS":
        return {"status": "RDF_UNRESOLVED", "reason": cn_point.get("reason"), "point": cn_point}
    r, raw_g = cn_point["r_A"], cn_point["raw_g"]
    peak_idx = int(cn_point["peak_index"])
    valley_A = float(cn_point["valley_A"])
    cn_series = np.asarray([
        _coordination_number(r, normalized_rdf(f.rdf_counts, 1, lock)[1], valley_A, lock)
        for f in frames
    ], dtype=float)
    g_peak_series = np.asarray([normalized_rdf(f.rdf_counts, 1, lock)[1][peak_idx] for f in frames], dtype=float)
    series["coordination_number"] = _finite_array(cn_series, "frame CN series")
    series["g_at_first_peak_bin"] = _finite_array(g_peak_series, "first-peak-bin g(r) series")
    taus = {name: integrated_autocorrelation_time(values) for name, values in series.items()}
    tau_max = max(taus.values())
    block_duration = max(0.5, 5.0 * tau_max)
    block_frames = int(math.ceil(block_duration / FRAME_INTERVAL_PS - 1e-12))
    effective_blocks = len(frames) // block_frames
    boot = lock["rdf"]["bootstrap"]
    if effective_blocks < int(boot["minimum_effective_blocks"]):
        return {
            "status": "RDF_UNCERTAINTY_INSUFFICIENT", "reason": "fewer than eight complete autocorrelation-sized blocks",
            "frame_count": len(frames), "tau_int_ps_by_observable": taus,
            "tau_int_max_ps": tau_max, "block_duration_ps": block_duration,
            "block_frames": block_frames, "complete_blocks": effective_blocks,
            "rdf_point": {k: v for k, v in cn_point.items() if k not in {"r_A", "raw_g", "smoothed_g"}},
        }
    reps = int(boot["replicates"])
    rng = np.random.default_rng(int(boot["seed"]))
    indices = _moving_block_indices(len(frames), block_frames, reps, rng)
    counts = np.asarray([f.rdf_counts for f in frames], dtype=np.int64)
    delta = np.full(reps, np.nan)
    cn_boot = np.full(reps, np.nan)
    g_boot = np.empty((reps, len(r)), dtype=float)
    detected = 0
    for b, idx in enumerate(indices):
        point = _rdf_point(counts[idx], lock)
        g_boot[b] = point["raw_g"] if point.get("status") == "PASS" else normalized_rdf(counts[idx].sum(axis=0), len(idx), lock)[1]
        if point.get("status") == "PASS":
            detected += 1
            delta[b] = point["peak_valley_delta_raw"]
            cn_boot[b] = point["coordination_number"]
    detect_fraction = detected / reps
    min_detect = 0.95
    delta_valid = delta[np.isfinite(delta)]
    cn_valid = cn_boot[np.isfinite(cn_boot)]
    if detect_fraction < min_detect or len(delta_valid) < int(math.ceil(min_detect * reps)):
        return {
            "status": "RDF_UNRESOLVED", "reason": "peak and valley detected in fewer than 95% of moving-block replicates",
            "bootstrap_detection_fraction": detect_fraction, "required_detection_fraction": min_detect,
        }
    delta_lower = float(np.quantile(delta_valid, 0.05))
    cn_ci = [float(np.quantile(cn_valid, 0.025)), float(np.quantile(cn_valid, 0.975))]
    windows = [_window_point(frames, bounds, lock) for bounds in final_windows_ps]
    cn_range = lock["rdf"]["coordination_number"]["accepted_range"]
    rdf_pass = bool(
        delta_lower > 0.0
        and cn_range[0] <= cn_point["coordination_number"] <= cn_range[1]
        and cn_ci[1] >= cn_range[0] and cn_ci[0] <= cn_range[1]
        and len(windows) == 2 and all(w.get("point_gate_pass") is True for w in windows)
    )
    smooth = cn_point["smoothed_g"]
    widths = peak_widths(smooth, [peak_idx], rel_height=0.5)[0]
    valley_idx = int(cn_point["valley_index"])
    higher, _ = find_peaks(smooth)
    higher = [int(i) for i in higher if r[i] > valley_A]
    band_low = np.quantile(g_boot, 0.025, axis=0).tolist()
    band_high = np.quantile(g_boot, 0.975, axis=0).tolist()
    return {
        "status": "PASS" if rdf_pass else "FAIL",
        "rdf_gate_pass": rdf_pass,
        "frame_count": len(frames),
        "radial_grid_A": r.tolist(), "pooled_raw_g": raw_g.tolist(), "pooled_smoothed_g": smooth.tolist(),
        "pooled_bootstrap_g_ci95_low": band_low, "pooled_bootstrap_g_ci95_high": band_high,
        "peak_A": cn_point["peak_A"], "valley_A": valley_A,
        "peak_height_raw": cn_point["peak_height_raw"], "valley_height_raw": cn_point["valley_height_raw"],
        "peak_valley_delta_raw": cn_point["peak_valley_delta_raw"],
        "peak_valley_delta_one_sided_95pct_lower": delta_lower,
        "peak_valley_bootstrap_detection_fraction": detect_fraction,
        "coordination_number": cn_point["coordination_number"], "coordination_number_ci95": cn_ci,
        "final_two_window_point_checks": windows,
        "tau_int_ps_by_observable": taus, "tau_int_max_ps": tau_max,
        "moving_block_duration_ps": block_duration, "moving_block_frames": block_frames,
        "complete_effective_blocks": effective_blocks, "bootstrap_replicates": reps,
        "first_peak_fwhm_A_report_only": float(widths[0] * (r[1] - r[0])),
        "higher_peak_positions_and_heights_report_only": [{"r_A": float(r[i]), "g_smoothed": float(smooth[i])} for i in higher],
        "valley_index": valley_idx,
    }


def _window_order_gate(frames: Sequence[FrameMetrics], bounds: Sequence[float], lock: dict = LOCK) -> dict:
    start, end = map(float, bounds)
    selected = [f for f in frames if start <= f.time_ps < end]
    expected = int(round((end - start) / FRAME_INTERVAL_PS))
    if len(selected) < expected - 1:
        return {"status": "QA_HOLD", "reason": "final 2 ps structural window is incomplete", "frame_count": len(selected)}
    q6_lock = lock["order"]["legacy_global_q6"]
    match = re.search(r"global_Q6 < ([0-9]+(?:\.[0-9]+)?)", q6_lock["joint_liquid_gate"]["per_frame_condition"])
    if not match:
        return {"status": "QA_HOLD", "reason": "unrecognized locked global-Q6 threshold"}
    q_threshold = float(match.group(1))
    q_values = [f.global_q6 for f in selected]
    if any(value is None or not math.isfinite(float(value)) for value in q_values):
        return {"status": "QA_HOLD", "reason": "global Q6 missing or non-finite"}
    q_fraction = float(np.mean([float(value) < q_threshold for value in q_values]))
    ptm_lock = lock["order"]["other_crystalline_order_gates"]
    ptm_limit = float(ptm_lock["crystalline_atom_fraction_max"])
    cluster_limit = int(ptm_lock["largest_connected_crystalline_cluster_atoms_max"])
    ptm_fraction = float(np.mean([
        f.ptm_crystalline_fraction is not None and f.ptm_largest_cluster_atoms is not None
        and f.ptm_crystalline_fraction <= ptm_limit and f.ptm_largest_cluster_atoms <= cluster_limit
        for f in selected
    ]))
    req = float(q6_lock["joint_liquid_gate"]["minimum_passing_frame_fraction"])
    req_ptm = float(ptm_lock["minimum_fraction_of_frames_meeting_limits_each_final_2ps_window"])
    passed = q_fraction >= req and ptm_fraction >= req_ptm
    return {
        "status": "PASS" if passed else "FAIL", "frame_count": len(selected),
        "global_q6_frame_pass_fraction": q_fraction, "global_q6_fraction_required": req,
        "mean_local_q6_report_only": float(np.mean([f.mean_local_q6 for f in selected])),
        "ptm_cluster_pass_fraction": ptm_fraction, "ptm_fraction_required": req_ptm,
        "q6_threshold": q_threshold, "ptm_crystalline_atom_fraction_max": ptm_limit,
        "largest_cluster_atoms_max": cluster_limit,
    }


def analyze_msd(unwrapped_positions_A: Sequence[Sequence[Sequence[float]]], block_frames: int,
                lock: dict = LOCK, replicates: int = BOOTSTRAP_REPLICATES) -> dict:
    x = _finite_array(unwrapped_positions_A, "unwrapped positions")
    n_atoms = int(lock["finite_cell"]["atoms"])
    if x.ndim != 3 or x.shape[1:] != (n_atoms, 3) or len(x) < 81:
        return {"status": "QA_HOLD", "reason": "MSD needs the complete 4 ps QA interval at 0.05 ps frame spacing"}
    dt = FRAME_INTERVAL_PS
    max_lag = len(x) - 1
    lags = np.arange(1, max_lag + 1)
    msd = np.asarray([
        np.mean(np.sum((x[k:] - x[:-k]) ** 2, axis=2)) for k in lags
    ], dtype=float)
    delta = float(msd[-1] - 0.0)
    # Strategy's frozen slope estimator: upper half of lags up to half of the
    # 4 ps QA interval, using multiple time origins.
    slope_max_frames = min(max_lag, int(round(2.0 / dt)))
    slope_start = max(1, slope_max_frames // 2)
    fit_lags = np.arange(slope_start, slope_max_frames + 1)
    if len(fit_lags) < 3:
        return {"status": "QA_HOLD", "reason": "insufficient lag range for MSD slope"}
    point_slope = float(np.polyfit(fit_lags * dt, msd[fit_lags - 1], 1)[0])
    increments = np.diff(x, axis=0)
    if len(increments) // block_frames < int(lock["rdf"]["bootstrap"]["minimum_effective_blocks"]):
        return {"status": "QA_HOLD", "reason": "MSD has fewer than eight complete moving blocks"}
    rng = np.random.default_rng(int(lock["rdf"]["bootstrap"]["seed"]))
    idx = _moving_block_indices(len(increments), block_frames, int(replicates), rng)
    slopes = np.empty(int(replicates), dtype=float)
    for b, inds in enumerate(idx):
        path = np.concatenate((np.zeros((1, n_atoms, 3)), np.cumsum(increments[inds], axis=0)), axis=0)
        b_msd = np.asarray([
            np.mean(np.sum((path[k:] - path[:-k]) ** 2, axis=2)) for k in fit_lags
        ], dtype=float)
        slopes[b] = float(np.polyfit(fit_lags * dt, b_msd, 1)[0])
    lower = float(np.quantile(slopes, 0.05))
    req_delta = float(lock["msd"]["pass"].split("delta MSD >= ")[1].split(" A^2")[0])
    passed = delta >= req_delta and lower > 0.0
    return {
        "status": "PASS" if passed else "FAIL", "msd_gate_pass": passed,
        "delta_msd_A2": delta, "delta_msd_lag_ps": float(max_lag * dt),
        "multiple_origin_slope_A2_per_ps": point_slope,
        "long_time_slope_one_sided_95pct_lower_A2_per_ps": lower,
        "msd_lags_ps": (lags * dt).tolist(), "msd_A2": msd.tolist(),
        "bootstrap_replicates": int(replicates), "moving_block_frames": int(block_frames),
    }


def _thermal_gate(first: Sequence[FrameMetrics], second: Sequence[FrameMetrics], target_K: float,
                  block_frames: int, thresholds: dict) -> dict:
    if len(first) < 2 or len(second) < 2:
        return {"status": "QA_HOLD", "reason": "thermal windows lack samples"}
    t1 = _finite_array([f.temperature_K for f in first], "first-window temperature")
    t2 = _finite_array([f.temperature_K for f in second], "second-window temperature")
    e1 = _finite_array([f.potential_energy_eV_per_atom for f in first], "first-window energy")
    e2 = _finite_array([f.potential_energy_eV_per_atom for f in second], "second-window energy")
    if len(t1) < block_frames or len(t2) < block_frames:
        return {"status": "QA_HOLD", "reason": "thermal uncertainty window is shorter than a correlation block"}
    rng = np.random.default_rng(int(LOCK["rdf"]["bootstrap"]["seed"]))
    reps = int(LOCK["rdf"]["bootstrap"]["replicates"])
    idx1 = _moving_block_indices(len(t1), block_frames, reps, rng)
    idx2 = _moving_block_indices(len(t2), block_frames, reps, rng)
    bmean1 = np.asarray([t1[i].mean() for i in idx1]); bmean2 = np.asarray([t2[i].mean() for i in idx2])
    ci1 = [float(np.quantile(bmean1, 0.025)), float(np.quantile(bmean1, 0.975))]
    ci2 = [float(np.quantile(bmean2, 0.025)), float(np.quantile(bmean2, 0.975))]
    all_t = np.concatenate((t1, t2))
    # The final two windows are adjacent and sampled on the locked 0.05 ps grid.
    slope = _slope(all_t, FRAME_INTERVAL_PS)
    means_tol = float(thresholds["temperature_mean_tolerance_percent"]) / 100.0
    adjacent_tol = float(thresholds["adjacent_window_mean_difference_tolerance_percent"]) / 100.0
    energy_delta = abs(float(e2.mean() - e1.mean()))
    passed = bool(
        abs(float(t1.mean()) - target_K) <= means_tol * target_K
        and abs(float(t2.mean()) - target_K) <= means_tol * target_K
        and abs(float(t2.mean() - t1.mean())) <= adjacent_tol * target_K
        and abs(slope) <= float(thresholds["absolute_temperature_trend_limit_K_per_ps"])
        and ci1[0] <= target_K <= ci1[1] and ci2[0] <= target_K <= ci2[1]
        and energy_delta <= float(thresholds["potential_energy_per_atom_adjacent_window_difference_limit_eV"])
    )
    return {
        "status": "PASS" if passed else "FAIL", "temperature_gate_pass": passed,
        "temperature_mean_first_K": float(t1.mean()), "temperature_mean_second_K": float(t2.mean()),
        "temperature_ci95_first_K": ci1, "temperature_ci95_second_K": ci2,
        "temperature_slope_K_per_ps": slope,
        "energy_mean_first_eV_per_atom": float(e1.mean()), "energy_mean_second_eV_per_atom": float(e2.mean()),
        "energy_adjacent_delta_eV_per_atom": energy_delta,
    }


def assess_phase(frames: Sequence[FrameMetrics], final_windows_ps: Sequence[Sequence[float]],
                 unwrapped_positions_A: Sequence[Sequence[Sequence[float]]], target_K: float,
                 q6_regression: dict, thermal_thresholds: dict, lock: dict = LOCK,
                 pressure_stress: dict | None = None) -> dict:
    """Complete V2 liquid gate; every missing/unresolved component yields HOLD."""
    if not frames:
        return {"gate_status": "PENDING", "reason": "no eligible phase frames yet"}
    if q6_regression.get("regression_pass") is not True:
        return {"gate_status": "QA_HOLD", "reason": "LEGACY_GLOBAL_Q6_REGRESSION_HOLD", "q6_regression": q6_regression}
    if any(f.ptm_status != "PASS" or not f.pbc_verified for f in frames):
        statuses = sorted({f.ptm_status for f in frames if f.ptm_status != "PASS" or not f.pbc_verified})
        return {"gate_status": "QA_HOLD", "reason": "PTM_UNAVAILABLE_OR_PERIODICITY_UNVERIFIED", "ptm_statuses": statuses}
    for bounds in final_windows_ps:
        expected = int(round((float(bounds[1]) - float(bounds[0])) / FRAME_INTERVAL_PS))
        selected = [f for f in frames if float(bounds[0]) <= f.time_ps < float(bounds[1])]
        if len(selected) < expected - 1:
            return {"gate_status": "PENDING", "reason": "FINAL_WINDOWS_INCOMPLETE", "window_frames": len(selected), "window_expected": expected}
    rdf = analyze_rdf_phase(frames, final_windows_ps, lock)
    if rdf.get("status") == "RDF_UNCERTAINTY_INSUFFICIENT":
        return {"gate_status": "PENDING", "reason": rdf.get("status"), "rdf": rdf}
    if rdf.get("status") not in {"PASS", "FAIL"}:
        return {"gate_status": "QA_HOLD", "reason": rdf.get("status"), "rdf": rdf}
    window_order = [_window_order_gate(frames, w, lock) for w in final_windows_ps]
    if any(w.get("status") == "QA_HOLD" for w in window_order):
        return {"gate_status": "QA_HOLD", "reason": "STRUCTURE_WINDOW_DATA_INSUFFICIENT", "rdf": rdf, "order_windows": window_order}
    block_frames = int(rdf["moving_block_frames"])
    ptm = np.asarray([f.ptm_crystalline_fraction for f in frames], dtype=float)
    bootstrap_idx = _moving_block_indices(
        len(frames), block_frames, int(lock["rdf"]["bootstrap"]["replicates"]),
        np.random.default_rng(int(lock["rdf"]["bootstrap"]["seed"])),
    )
    local_q6_slope_ci = _bootstrap_ci(np.asarray([f.mean_local_q6 for f in frames], dtype=float), bootstrap_idx, "slope")
    ptm_slope_ci = _bootstrap_ci(ptm, bootstrap_idx, "slope")
    ordering_trend = ptm_slope_ci[0] > 0.0
    thermal = _thermal_gate(
        [f for f in frames if final_windows_ps[0][0] <= f.time_ps < final_windows_ps[0][1]],
        [f for f in frames if final_windows_ps[1][0] <= f.time_ps < final_windows_ps[1][1]],
        target_K, block_frames, thermal_thresholds,
    )
    if thermal.get("status") == "QA_HOLD":
        if thermal["reason"] in {"thermal windows lack samples", "thermal uncertainty window is shorter than a correlation block"}:
            return {"gate_status": "PENDING", "reason": thermal["reason"], "rdf": rdf, "order_windows": window_order}
        return {"gate_status": "QA_HOLD", "reason": thermal["reason"], "rdf": rdf, "order_windows": window_order}
    msd = analyze_msd(unwrapped_positions_A, block_frames, lock)
    if msd.get("status") == "QA_HOLD":
        if "needs the complete 4 ps QA interval" in msd["reason"] or "fewer than eight complete moving blocks" in msd["reason"]:
            return {"gate_status": "PENDING", "reason": msd["reason"], "rdf": rdf, "order_windows": window_order, "thermal": thermal}
        return {"gate_status": "QA_HOLD", "reason": msd["reason"], "rdf": rdf, "order_windows": window_order, "thermal": thermal}
    if not isinstance(pressure_stress, dict):
        stress_gate = pressure_stress_gate(None, None, lock)
    elif pressure_stress.get("input_error"):
        stress_gate = {"status": "QA_HOLD", "pressure_stress_gate_pass": False,
                       "reason": str(pressure_stress["input_error"])}
    else:
        stress_gate = pressure_stress_gate(
            pressure_stress.get("times_ps"), pressure_stress.get("stress_eV_A3"), lock,
        )
    if stress_gate["status"] == "PENDING":
        return {
            "gate_status": "PENDING", "reason": "PRESSURE_STRESS_DATA_INSUFFICIENT",
            "rdf": rdf, "order_windows": window_order, "thermal": thermal, "msd": msd,
            "pressure_stress": stress_gate,
        }
    if stress_gate["status"] == "QA_HOLD":
        return {
            "gate_status": "QA_HOLD", "reason": "PRESSURE_STRESS_QA_HOLD",
            "rdf": rdf, "order_windows": window_order, "thermal": thermal, "msd": msd,
            "pressure_stress": stress_gate,
        }
    order_pass = all(w.get("status") == "PASS" for w in window_order) and not ordering_trend
    passed = bool(
        rdf.get("rdf_gate_pass") and order_pass and thermal.get("temperature_gate_pass")
        and msd.get("msd_gate_pass") and stress_gate.get("pressure_stress_gate_pass")
    )
    if q6_regression.get("regression_pass") is not True:
        return {
            "gate_status": "QA_HOLD", "reason": "LEGACY_GLOBAL_Q6_REGRESSION_HOLD",
            "q6_regression": q6_regression, "rdf": rdf, "order_windows": window_order,
            "ordering_trend": ordering_trend, "mean_local_q6_slope_ci95_report_only": local_q6_slope_ci,
            "ptm_fraction_slope_ci95": ptm_slope_ci, "thermal": thermal, "msd": msd,
            "pressure_stress": stress_gate,
        }
    return {
        "gate_status": "PASS" if passed else "FAIL_WINDOW",
        "reason": (
            "PRESSURE_STRESS_DRIFT" if stress_gate["status"] == "FAIL_DRIFT" else
            "ORDERING_TREND" if ordering_trend else None
        ),
        "rdf": rdf, "order_windows": window_order,
        "ordering_trend": ordering_trend,
        "mean_local_q6_slope_ci95_report_only": local_q6_slope_ci,
        "ptm_fraction_slope_ci95": ptm_slope_ci,
        "thermal": thermal, "msd": msd,
        "pressure_stress": stress_gate,
        "ptm_implementation": "OVITO PolyhedralTemplateMatchingModifier",
        "ptm_version": frames[-1].ptm_version, "ptm_rmsd_cutoff": 0.10,
    }


def candidate_snapshot_spacing_steps(frames: Sequence[FrameMetrics], thermal_values: Sequence[float] | None = None,
                                     energy_values: Sequence[float] | None = None, lock: dict = LOCK) -> dict:
    if len(frames) < 8 or any(f.ptm_status != "PASS" for f in frames):
        raise CriteriaError("candidate spacing requires at least eight complete PTM-qualified frames")
    # Use the same trajectory-frame cadence (0.05 ps) for every observable.
    # Thermal/energy logs may be more frequent and must not be mixed onto this
    # grid.  The pooled first-valley location supplies a consistent per-frame CN.
    pooled = _rdf_point(np.asarray([f.rdf_counts for f in frames], dtype=np.int64), lock)
    if pooled.get("status") != "PASS":
        raise CriteriaError("candidate spacing requires a resolved pooled RDF first valley")
    r = pooled["r_A"]
    valley = float(pooled["valley_A"])
    cn_values = []
    for frame in frames:
        _r, g = normalized_rdf(frame.rdf_counts, 1, lock)
        cn_values.append(_coordination_number(_r, g, valley, lock))
    fields = {
        "temperature": [f.temperature_K for f in frames],
        "potential_energy_per_atom": [f.potential_energy_eV_per_atom for f in frames],
        "mean_local_q6": [f.mean_local_q6 for f in frames],
        "ptm_crystalline_fraction": [f.ptm_crystalline_fraction for f in frames],
        "coordination_number": cn_values,
    }
    taus = {name: integrated_autocorrelation_time(_finite_array(values, name)) for name, values in fields.items()}
    tau_max = max(taus.values())
    spacing_ps = max(0.5, 5.0 * tau_max)
    spacing_steps = int(math.ceil(spacing_ps / 0.05) * STEP_INTERVAL)
    return {
        "tau_int_ps_by_observable": taus, "tau_int_max_ps": tau_max,
        "minimum_snapshot_spacing_ps": spacing_ps,
        "snapshot_spacing_steps": spacing_steps,
        "snapshot_spacing_ps_rounded_to_frame_grid": spacing_steps * 0.0005,
        "source_family_policy": "all snapshots remain in this single task family",
    }
