#!/usr/bin/env python3
"""Map one QSGS seed to 50/100/150³ and measure mixed-cell/source effects."""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
PHASE = PROJECT / "phase_change_validation"
sys.path.insert(0, str(PHASE))
sys.path.insert(0, str(HERE))

from diagnostics import mapping_diagnostics, write_csv
from material_model import CONFIG
from run_3d_phase_change import load_seed_ag


def main():
    seed = load_seed_ag()
    cfg = CONFIG["numerical_heating_test"]
    edge = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    rows = []
    for n in (50, 100, 150):
        rows.extend(mapping_diagnostics(
            seed, n, edge, float(cfg["source_rate_W_m3_at_10um"]),
            float(cfg["source_region_radius_fraction_of_edge"]),
        ))
    output = HERE / "results" / "QSGS_MAPPING_CONNECTIVITY.csv"
    write_csv(output, rows)
    print(f"wrote {len(rows)} rows to {output}")


if __name__ == "__main__":
    main()
