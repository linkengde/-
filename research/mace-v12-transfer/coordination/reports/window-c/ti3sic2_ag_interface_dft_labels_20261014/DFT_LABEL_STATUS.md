# Stage C-D1 DFT label status

```text
DFT_INPUT_READY=YES
DFT_LABEL_GENERATION_RUNNING=NO
DFT_LABEL_READY=NO
STRUCTURES_AUDITED=15/15
SCF_STARTED=0
DFT_LABELS_CREATED=0
```

The 15 geometries and configs are bound to input package commit `1fec365ecedb8555ae2c857509d1018953639b86`. All 15 structure/config SHA checks passed. No new DFT result exists; blank energy/force/stress fields are missing labels, not zero-valued data.

## Stop reason

No GPAW job was launched because the current execution container fails the safe resource gate. At the recorded snapshot it has a four-CPU quota, 34,359,738,368-byte (32 GiB) cgroup memory limit and 28,322,897,920 bytes in use, leaving 6,036,840,448 bytes under the cgroup cap. Its recorded cgroup high-water events include 15 `max` events, although current `oom` and `oom_kill` counters are zero. `/workspace` has 22,072,643,584 bytes free. The previous E7 resource assessment explicitly rated this environment `RESOURCE_AVAILABLE=NO` and instructed not to launch even the 48-atom calculations here without a fresh task-specific reserve/peak-memory basis.

The 217-atom interface scan also has no demonstrated 2×2×1 k-point convergence. No denser-mesh result or error budget is present. The nine interface calculations therefore have both resource and method-convergence holds. Precision parameters were not weakened.

No TSCBOP file was modified, no fitting was run, and no MD or TTM was launched.
