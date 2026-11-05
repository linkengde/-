#!/usr/bin/env python3
"""Synthetic regression test for the P1G VOF upper-envelope surface height."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


MODULE = Path(__file__).with_name("p1g_postprocess.py")
spec = importlib.util.spec_from_file_location("p1g_postprocess", MODULE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

fields = np.zeros((20, 2, 2, 4), dtype=np.float32)
fields[0, 0, 0, [0, 2]] = 1.0  # internal gap: highest occupied cell must win
fields[0, 0, 1, 3] = 1.0
fields[5] = np.arange(4, dtype=np.float32)[None, None, :] + 300.0
meta = {"dx_m": 1.0e-9, "origin_m": (0.0, 0.0, 5.0e-9)}

_, _, height, surface_temperature = module.derive(fields, meta)
assert np.isclose(height[0, 0], 7.5e-9)
assert np.isclose(height[0, 1], 8.5e-9)
assert np.isnan(height[1, 1])
assert surface_temperature[0, 0] == 302.0
assert surface_temperature[0, 1] == 303.0
print("P1G upper-envelope surface regression passed")
