#!/usr/bin/env python3
"""Export a v11 wrapper whose serialized MACE config explicitly selects CPU."""
from pathlib import Path
import sys
import torch

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = Path(__file__).resolve().parent
MODEL = ROOT / "checkpoints/MACE_periodic_v11_energyfocus_run-44.model"
OUT = OUTDIR / "mace_v11_lammps_mliap_forcecpu.pt"
INTEGRATION = Path("/workspace/.local/lammps-ttm/run_1ps_diagnostic_20261001/potential_screening/lammps_mace_integration")
sys.path.insert(0, str(INTEGRATION))
from mace.calculators.lammps_mliap_mace import LAMMPS_MLIAP_MACE

if not MODEL.exists():
    raise SystemExit(f"v11 model not available: {MODEL}")
if OUT.exists():
    raise SystemExit(f"Refusing to overwrite {OUT}")

model = torch.load(MODEL, map_location="cpu", weights_only=False).double().cpu()
model.lammps_mliap = True
wrapper = LAMMPS_MLIAP_MACE(model, head="Default")
wrapper.config.force_cpu = True
torch.save(wrapper, OUT)
print(f"Saved {OUT}")
print(f"Elements: {wrapper.element_types}")
print(f"force_cpu: {wrapper.config.force_cpu}")
print(f"Cutoff: {float(model.r_max):.6f} A")
