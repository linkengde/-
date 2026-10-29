# Ag70 1× large-pool cooling energy audit

## Finding

For the Ag 70 vol% 1× Case A model, the 19.25 ps heat-off state cannot reach an all-solid equilibrium while the existing assumptions remain in place. The heat-off enthalpy is **2.908798024218367×10⁻¹⁴ J**, above the model's all-Ag-solid-at-(T_m) enthalpy **(E_s=2.611689713605194×10^{-14}) J**. At a hypothetical uniform-(T_m) equilibrium, energy balance leaves an estimated **26.9305% of Ag liquid**.

Cooling was safely stopped at the latest complete and independently restored checkpoint, **19.49 ps**. The thermal history was not restarted or altered, and no data were removed. The result is an energy-balance conclusion within the candidate-property Case A model, not a validated prediction for the physical composite.

## Scope and energy definition

This audit uses the existing Ag70 QSGS geometry (actual grid Ag fraction 69.9996 vol%), the 1× single-temperature enthalpy solver, Case A, zero Ti₃SiC₂ heat capacity and conductivity, and insulated external faces. No 2×, TTM, or new heating run was started for this audit.

All energies below are enthalpy relative to the solver's 300 K reference. (E_s) is the sum of the temperature-dependent Ag solid sensible enthalpy evaluated at (T_m), with no latent heat. (E_l=E_s+E_{latent}) adds the latent heat for all Ag. These values are integrated over the actual Ag voxels and their native cell volume.

The candidate inputs used by the solver are: Ag density 10,490 kg m⁻³; (T_m=1235.08) K; latent heat 104,729.66 J kg⁻¹; solid (C_p) anchors 235.47 J kg⁻¹ K⁻¹ at 300 K and 294.80 J kg⁻¹ K⁻¹ at (T_m); liquid (C_p=310.29) J kg⁻¹ K⁻¹ at its single anchor. The solid (C_p) is linearly interpolated between its anchors, and the liquid value is held constant by the current model. These candidate properties are not fully source-verified.

The Ag volume represented by the grid is 1.0042104138186939×10⁻²³ m³, corresponding to 1.0534167240958099×10⁻¹⁹ kg with the candidate density.

## Energy comparison

| Quantity | Value (J) | Meaning |
| --- | ---: | --- |
| (E_{off}) | 2.908798024218367×10⁻¹⁴ | Total enthalpy in the 19.25 ps checkpoint when heating ended |
| (E_s) | 2.611689713605194×10⁻¹⁴ | All Ag solid at (T_m), just before latent melting |
| (E_l) | 3.7149294758437745×10⁻¹⁴ | All Ag liquid at (T_m) |
| (E_{off}-E_s) | 2.9710831061317302×10⁻¹⁵ | Excess over the all-solid-at-(T_m) enthalpy |
| (E_l-E_{off}) | 8.061314516254074×10⁻¹⁵ | Energy still needed to melt all Ag at (T_m) |

Thus (E_s<E_{off}<E_l). The energy-based uniform-(T_m) estimate is

\[
f_{eq}=\frac{E_{off}-E_s}{E_l-E_s}=0.2693053\quad(26.9305\%).
\]

This is an idealized energy-balance estimate, not evidence that the finite-time field has reached equilibrium. It nevertheless proves that the present insulated Case A enthalpy model cannot convert all Ag to solid while conserving the measured (E_{off}): its all-solid enthalpy ceiling is (E_s<E_{off}).

## Cooling checks

The source-off driver calls `model.run(source_on=False)`. The cumulative source-energy ledger remains constant at 2.9087980241873447×10⁻¹⁴ J through the cooling records. At the latest checkpoint, total enthalpy is 2.908798024218365×10⁻¹⁴ J; the residual is 3.102053514818689×10⁻²⁵ J, or 1.0664×10⁻¹¹ relative to the injected energy.

| Time (ps) | Peak T (K) | Ag-weighted mean T (K) | Ag liquid fraction | Largest connected cluster (nm³) |
| ---: | ---: | ---: | ---: | ---: |
| 18.75 heat-off | 2300.06 | 1248.31 | 22.7977% | 2285.05 |
| 19.25 | 1394.72 | 1237.54 | 26.0776% | 2613.40 |
| 19.27 | 1386.81 | 1237.29 | 26.1553% | 2621.29 |
| 19.49 latest complete | 1324.43 | 1235.22 | 26.7933% | 2685.93 |

The 19.25, 19.27, and 19.49 ps checkpoints all passed `cool_ag70_from_checkpoint.py --restore-only`. The maximum connected liquid cluster remains one face-connected component at these samples. The liquid fraction rises as the peak temperature falls because stored heat redistributes from the hottest, superheated voxels into neighboring Ag. More voxels enter the enthalpy model's latent-heat plateau near (T_m). The Ag-weighted mean temperature approaches (T_m), while the hottest voxel remains above it because the field is still spatially nonuniform at 19.49 ps.

The heat-off step peaked at 2300.05995 K, about 0.06 K above the configured 2300 K numerical guard due to the final integration step. Vaporization is absent, so this high-temperature cooling branch remains exploratory.

## Decision on further cooling

No further long-duration Case A cooling is needed to answer whether this heat-off state can become fully solid: energy conservation rules that out under the current model. The run was stopped after preserving and validating the 19.49 ps checkpoint. It was not forced to zero liquid fraction.

For a physically meaningful full-freezing path, the model needs real Ti₃SiC₂ heat capacity and conductivity, an Ag/Ti₃SiC₂ interface treatment, and a physically justified heat-rejection path such as conduction into the contact/electrode bulk. Adding ceramic storage can absorb energy internally; because that alone does not remove total system energy, a finite external thermal boundary is also needed when the composite's conserved energy exceeds its all-solid capacity. Property sources and boundary parameters must be validated before quantitative solidification times are claimed.

## Reproducibility and files

- Analysis script: [audit_ag70_1x_large_pool_cooling.py](audit_ag70_1x_large_pool_cooling.py)
- Checkpoint metrics: [cooling_energy_checkpoint_metrics.csv](cooling_energy_checkpoint_metrics.csv)
- Cooling trace: [cooling_trace_18p75ps_to_latest.csv](cooling_trace_18p75ps_to_latest.csv)
- Cooling plot: [cooling_trace_18p75ps_to_latest.png](cooling_trace_18p75ps_to_latest.png)
- Machine-readable results: [cooling_energy_audit_results.json](cooling_energy_audit_results.json)
- Latest validated checkpoint: [caseA_ag70_1x_cool_19p49ps.npz](../ag70vol_cooling_extend_from19p25_20261029/checkpoints/caseA_ag70_1x_cool_19p49ps.npz)

Recreate the tables, JSON, and plot from the repository root with:

```bash
python models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/ag70vol_cooling_energy_audit_20261029/audit_ag70_1x_large_pool_cooling.py
```
