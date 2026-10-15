# Pretrained potential selection decision

## Decision

Neither pretrained model is approved as a replacement production potential. MatterSim-v1.0.0-5M is the more promising frozen comparison for interface forces on this one 135-atom geometry, but it is not physically consistent across the checks: its composition-subtracted interaction energy is 7.269 eV less binding than DFT native-total convention (6.419 eV relative to the distinct force-consistent-free convention), and its net normal interface force has the opposite sign and much larger magnitude. MatterSim's own model card cautions against quantitative surface/interface use without fine-tuning. SevenNet has no overall improvement over the best existing MACE full-interface force baseline and retains substantial Si/Ti and cross-force errors.

MatterSim full-interface force RMSE `0.18489 eV/Å` is 22.6% lower than the best existing MACE full-interface result in the table (V17 `0.23885`) and 24.0% lower than MACE-MP (`0.24326`). It is also lower on the 114-atom TSC-only component. However, the 48-atom bulk control is `0.04562` vs MACE-MP `0.03817`, and the interface energy/normal-force failures invalidate a general “more reasonable” conclusion. SevenNet full-interface RMSE `0.26282 eV/Å` is worse than MACE-MP and V17.

## Per-model gate table

| Model | MODEL_LOAD_PASS | DFT_STATIC_RESULT | TSC_BACKBONE_RESULT | INTERFACE_FORCE_RESULT | COMPUTATIONAL_FEASIBILITY | IMPROVEMENT_OVER_MACE | SMALL_MD_READINESS | PRODUCTION_READY |
|---|---|---|---|---|---|---|---|---|
| SevenNet-MF-0 (PBE) | YES | YES (finite on 48 and 135/component cases) | CONDITIONAL: bulk RMSE 0.05676; TSC-slab 0.40616, no better than V17 | NO overall: full 0.26282 and cross 0.41841; normal-force sign/magnitude disagree | YES for CPU inference at screened size | NO | NO | NO |
| MatterSim-v1.0.0-5M | YES | YES (finite on 48 and 135/component cases) | CONDITIONAL: slab RMSE 0.20586 improves; 48 bulk 0.04562 slightly worse than MACE-MP | CONDITIONAL/NO as overall interface evidence: full 0.18489/cross 0.18459 improve force RMSE, but interaction energy and normal response fail | YES for CPU inference at screened size | PARTIAL only; not an overall win | NO | NO |

## Data independence and limits

The pre-existing project evaluation inventory records MACE V12 as 23 train / 2 validation / 3 test structures; V14 evaluation uses the run-45 checkpoint and its earlier gate is SCREEN_FAIL; V17 has 43 train / 3 reused development frames and no test file; V18 has 30 train / 3 development frames and no test file. The inherited exact-geometry scan found no exact match for the 135-atom full geometry or its Ag/TSC subsets in the available V12/V14/V17/V18 partitions. A parent-family join across all historically scored development data is incomplete, so family-level independence is UNVERIFIED; this is not a blind test. V18's 48-atom backbone failure remains material because its training set contained no pure Ti₃SiC₂ examples.

The 135-atom interface, Ag-only and TSC-only DFT data are same-geometry references from one family, and the 48-atom bulk is one additional control. None is a broad independent blind test. Exact/family overlap with SevenNet's Materials Project data and MatterSim's training corpus has not been established. The three gap structures from the prior task were not rerun in this screen and would also be one parent family without DFT labels. Thus no model may be declared production-ready or ready for MD from this task.

Recommended next step is to retain MACE-MP/V17 as project comparison baselines and MatterSim-v1.0.0-5M as an additional static diagnostic only. If a model is to be considered for small MD, first obtain independent DFT interface structures covering distinct terminations/registries and distances, inspect matching force and energy increments, then obtain explicit authorization for MD. This task performed no new DFT, training/fine-tuning, relaxation, MD, or TTM.
