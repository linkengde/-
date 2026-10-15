#!/usr/bin/env python3
"""Export the final v11 MACE model to LAMMPS ML-IAP's Python wrapper format."""
from pathlib import Path
import sys
import torch

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = Path(__file__).resolve().parent
MODEL = ROOT / "checkpoints/MACE_periodic_v11_energyfocus_run-44.model"
OUT = OUTDIR / "mace_v11_lammps_mliap.pt"
INTEGRATION = Path("/workspace/.local/lammps-ttm/run_1ps_diagnostic_20261001/potential_screening/lammps_mace_integration")
sys.path.insert(0, str(INTEGRATION))
from mace.calculators.lammps_mliap_mace import LAMMPS_MLIAP_MACE

if not MODEL.exists():
    raise SystemExit(f"v11 model not yet available: {MODEL}")
if OUT.exists():
    raise SystemExit(f"Refusing to overwrite {OUT}")

model = torch.load(MODEL, map_location="cpu", weights_only=False).double().cpu()
model.lammps_mliap = True
wrapper = LAMMPS_MLIAP_MACE(model, head="Default")
torch.save(wrapper, OUT)
print(f"Saved {OUT}")
print(f"Elements: {wrapper.element_types}")
print(f"Cutoff: {float(model.r_max):.6f} A")
