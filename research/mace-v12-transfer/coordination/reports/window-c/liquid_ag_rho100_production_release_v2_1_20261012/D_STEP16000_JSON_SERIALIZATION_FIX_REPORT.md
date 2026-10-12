# Step 16000 JSON serialization recovery

The runner completed integration through step 16000 and saved a valid checkpoint, then stopped while writing the QA monitor. The monitor's RDF `final_two_window_point_checks` contained NumPy arrays (`r_A`, `raw_g`, `smoothed_g`), which `json.dump` cannot serialize. The corrected runner converts NumPy arrays/scalars to JSON-native lists/scalars while retaining `allow_nan=False`, so invalid non-finite results still fail closed.

The fix changes only serialization. It does not alter the QA criteria, scientific locks, EAM potential, cell, structure, velocities, RNG, timestep, thermostat, or temperature path. All 21 runner tests pass, including a new regression test for NumPy diagnostics and NaN rejection. The actual step 16000 monitor was rebuilt from the existing trajectory/logs and now serializes successfully as `PENDING` (`MSD has fewer than eight complete moving blocks`).

The immutable source checkpoint remains step 16000, SHA-256 `dfd90239e32e543d638d080f633cdf6b4b52392e0642c75623b83db91a578ddf`. Its positions, momenta, and all three RNG streams restore successfully. The incident snapshot is at `fad4502a713822857de7f15b0cdf9c5616710d04`; no MD ran after step 16000.

The new bounded continuation config SHA-256 is `c8cf22354e6ac1be40a244a1d4a02e7fad5b5176f2d5d313e874261a54786de6` and the direct user override SHA-256 is `d3c96f47c7acb8f1fce28cd5676c63b99c3e4afac4bdfdd58f2385e0c028a2ca`. It resumes from step 16000, caps total execution at step 60000, and retains 1000-step/3600-second segment limits. `execution_authorized` remains false in the preparation config; the direct user instruction is recorded separately.
