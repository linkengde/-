# CHECKPOINT_VALIDATION

Each checkpoint atomically writes double-precision temperature and volumetric enthalpy arrays to compressed NPZ. Its metadata stores physical time, next dt, source and boundary energy totals, initial stored enthalpy, step count, phase/phase time, last accepted-step profile, Newton profiles, material/solver settings, and geometry/source hashes. Static Ag fraction and source-rate fields are saved beside the checkpoint. Accepted-step JSONL is flushed and `fsync`ed before checkpoint replacement; failed trial reasons are separately flushed.

Restore verifies shape, float64 dtype, phase-mask SHA-256, source-field SHA-256, Ag melt settings, latent heat, Ti₃SiC₂ temperature guard, and solver tolerances/options before replacing state.

## Validation

`test_checkpoint_restart_matches_continuous_advance` advanced a deterministic 9×8×7 mixed-phase case, saved and restored after one step, then advanced both restored and uninterrupted simulations. Temperature and enthalpy arrays matched bit-for-bit; time, source energy, and energy accounting also matched. This verifies the restart path, including cumulative enthalpy/energy state.

The production runs also left checkpoints at each accepted step. The completed 100³ safe run is restartable at its final state, and the 150³ short test is restartable at its 16-step limit. The 50³ safe run is restartable at the cooled endpoint. After deduplicating identical static geometry/source arrays into `results/checkpoints/shared_inputs/`, all three representative saved states were reconstructed from the generated inputs and restored successfully. The N=50 restore was checked with BLAS threads fixed at 1 because conservative remapping involves floating-point reductions; the benchmark runner now enforces that setting before importing NumPy.

## Resume example

Use the exact checkpoint directory and original flags, then add `--resume`. A mismatch in geometry, source, thermophysical candidate settings, solver tolerances, resolution, policy, thread count, preconditioner, or time-step factor is rejected. Do not point two concurrent jobs at the same checkpoint directory.
