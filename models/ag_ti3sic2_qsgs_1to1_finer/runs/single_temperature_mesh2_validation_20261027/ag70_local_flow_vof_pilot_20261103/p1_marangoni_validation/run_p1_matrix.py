#!/usr/bin/env python3
"""Compile and run the isolated dimensionless P1 planar Marangoni matrix."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "p1_planar_marangoni.c"
CASES = HERE / "cases"
QCC = Path(os.environ.get("BASILISK_QCC", "/tmp/basilisk-build/src/qcc"))
BASILISK = Path(os.environ.get("BASILISK_SRC", str(QCC.parent)))
RESULTS = HERE / "results"
T_END = float(os.environ.get("P1_T_END", "5"))


def run(cmd: list[str], *, cwd: Path = HERE, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=cwd, env=env, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, check=False)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    CASES.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["BASILISK"] = str(BASILISK)
    env["PATH"] = f"{BASILISK}:{env.get('PATH', '')}"
    build = run([str(QCC), "-grid=quadtree", "-O2", SOURCE.name,
                 "-o", "p1_planar_marangoni", "-lm"], env=env)
    (RESULTS / "compile.log").write_text(build.stdout)
    if build.returncode:
        raise SystemExit(f"qcc failed ({build.returncode}); see results/compile.log")

    specs = [(lev, mode, 0.5) for lev in (4, 5, 6)
             for mode in ("zero", "neg", "pos", "agred", "agoxygen")]
    specs.append((5, "neg", 0.25))
    results: list[dict[str, str]] = []
    samples: list[dict[str, str]] = []
    for level, mode, dtf in specs:
        tag = f"Hdx{1 << level}_{mode}_dt{dtf:g}"
        proc = run([str(HERE / "p1_planar_marangoni"), str(level), mode,
                    f"{dtf:g}", f"{T_END:g}"])
        (CASES / f"{tag}.log").write_text(proc.stdout)
        if proc.returncode:
            raise SystemExit(f"case {tag} failed ({proc.returncode}); see cases/{tag}.log")
        for line in proc.stdout.splitlines():
            m = re.match(r"^(RESULT|SAMPLE) (.*)$", line)
            if not m:
                continue
            fields = dict(re.findall(r"([A-Za-z0-9_]+)=([^ ]+)", m.group(2)))
            if m.group(1) == "RESULT":
                fields["case"] = tag
                results.append(fields)
            else:
                fields["case"] = tag
                samples.append(fields)
    write_csv(RESULTS / "p1_raw_results.csv", results)
    write_csv(RESULTS / "p1_raw_samples.csv", samples)

    velocity = []
    mass = []
    stress = []
    timestep = []
    for row in results:
        level = row["level"]
        mode = row["mode"]
        velocity.append({k: row.get(k, "") for k in (
            "case", "level", "H_over_dx", "mode", "dt_factor", "t", "tau",
            "d_sigma_dT", "grad_T", "u_interface", "u_exact", "velocity_error_pct",
            "u_rms_interface", "umax", "mass_drift_pct", "cells")})
        mass.append({k: row.get(k, "") for k in (
            "case", "level", "H_over_dx", "mode", "dt_factor", "t",
            "mass_drift_pct", "u_interface", "umax")})
        stress.append({k: row.get(k, "") for k in (
            "case", "level", "H_over_dx", "mode", "dt_factor", "tau",
            "force_per_length", "force_expected", "viscous_stress_jump",
            "stress_balance_error_pct", "shear_top", "shear_bottom", "mgp_i", "mgp_resa")})
        if row["case"] in ("Hdx32_neg_dt0.5", "Hdx32_neg_dt0.25"):
            timestep.append({k: row.get(k, "") for k in (
                "case", "dt_factor", "t", "u_interface", "u_exact", "velocity_error_pct",
                "mass_drift_pct", "stress_balance_error_pct", "mgp_i", "mgp_resa")})
    write_csv(RESULTS / "p1_velocity_convergence.csv", velocity)
    write_csv(RESULTS / "p1_mass_conservation.csv", mass)
    write_csv(RESULTS / "p1_surface_stress_balance.csv", stress)
    write_csv(RESULTS / "p1_timestep_sensitivity.csv", timestep)

    metadata = {
        "solver": "Basilisk C qcc, 2-D quadtree uniform Cartesian grid, centered Stokes + VOF + tension.h automatic curvature",
        "source_file": SOURCE.name,
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "qcc_path": str(QCC),
        "qcc_sha256": hashlib.sha256(QCC.read_bytes()).hexdigest(),
        "basilisk_header_sha256": {
            name: hashlib.sha256((BASILISK / name).read_bytes()).hexdigest()
            for name in ("tension.h", "iforce.h", "curvature.h", "two-phase.h")
        },
        "basilisk_release_or_darcs_revision": "not recorded in the installed build; executable and key-header hashes are recorded",
        "dimensionless": True,
        "ag_properties_or_geometry_used": False,
        "interface": "flat liquid layer below y=0, gas above y=0; periodic x; no-slip top and bottom",
        "rho_liquid": 1.0, "rho_gas": 1.0,
        "mu_liquid": 1.0, "mu_gas": 1.0,
        "sigma_constant": 0.01,
        "temperature_gradient": 0.01,
        "d_sigma_dT": {"negative": -0.5, "positive": 0.5, "zero": 0.0},
        "silver_literature_slope_sensitivity": {
            "normalization": "dimensionless beta*=0.5*(beta/beta_ref), beta_ref=+0.546 mN/(m K); grad_T*=0.01; tau*=beta* grad_T*",
            "agred": {"beta_mN_per_m_K": -0.298, "atmosphere": "Ar/He + 5% H2; 2.6 ppm H2O; apparent oxygen partial pressure varies with temperature"},
            "agoxygen": {"beta_mN_per_m_K": 0.546, "oxygen_partial_pressure_Pa": 100.0},
            "source": "Morohoshi et al., J. Japan Inst. Metals 72 (2008), 708-713, doi:10.2320/jinstmet.72.708",
            "meaning": "sign and flow-direction sensitivity only; not an Ag dimensional-velocity prediction",
        },
        "surface_shear": "tau=d_sigma_dT*grad_T; imposed force density tau*|grad(f)| e_x",
        "analytical_interface_velocity": "U=tau*H/(mu_liquid+mu_gas) for equal layer heights H; here H=mu_liquid=mu_gas=1, so U=tau/2",
        "H_over_dx": [16, 32, 64],
        "timestep_factors": [0.5, 0.25],
        "final_time": T_END,
        "official_basilisk_marangoni_case": "https://basilisk.fr/src/test/marangoni.c; not present in the installed source mirror, so no code-to-code reproduction is claimed",
        "surface_force_closure": "run_p1_closure.py integrates automatic-curvature normal and analytical tangential forces over PLIC facets of a static circle",
    }
    (RESULTS / "p1_solver_configuration.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"completed {len(results)} runs; results in {RESULTS}")


if __name__ == "__main__":
    main()
