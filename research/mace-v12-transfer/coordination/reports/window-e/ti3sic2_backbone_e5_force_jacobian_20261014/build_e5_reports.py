#!/usr/bin/env python3
"""Summarize verified Stage E5 run-0 force-response data; no simulation calls."""
import csv, json, math
from pathlib import Path
import numpy as np

ROOT=Path('/workspace/-')
OUT=ROOT/'research/mace-v12-transfer/coordination/reports/window-e/ti3sic2_backbone_e5_force_jacobian_20261014'
E4=ROOT/'research/mace-v12-transfer/coordination/reports/window-e/ti3sic2_backbone_e4_force_jacobian_20261014'
E3=ROOT/'research/mace-v12-transfer/coordination/reports/window-e/ti3sic2_backbone_e3_fit_design_20261014'
POT=ROOT/'research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/TSCBOP.tersoff'

def read_csv(path):
    with path.open(newline='') as f: return list(csv.DictReader(f))

def fmt(x, n=4):
    return f'{x:.{n}g}' if isinstance(x,(int,float,np.number)) and math.isfinite(float(x)) else 'n/a'

def main():
    status=json.loads((OUT/'run0_status.json').read_text())
    analysis=json.loads((OUT/'jacobian_analysis.json').read_text())
    matrix_rows=read_csv(OUT/'PARAMETER_RESPONSE_MATRIX.csv')
    cols=list(matrix_rows[0])[4:]
    J=np.array([[float(r[c]) for c in cols] for r in matrix_rows],dtype=float)
    cases=np.array([r['case'] for r in matrix_rows])
    fam=np.array(['bulk_parent' if c.startswith('bulk_') else 'derived_slab' for c in cases])
    fam_counts={k:int(np.sum(fam==k)) for k in sorted(set(fam))}
    weights=np.array([1/math.sqrt(fam_counts[x]) for x in fam])
    Jfb=J*weights[:,None]
    norms=np.linalg.norm(Jfb,axis=0); active=norms>max(norms)*1e-12
    Jfbu=Jfb[:,active]/norms[active]
    sfb=np.linalg.svd(Jfbu,compute_uv=False)
    analysis['family_balanced_unit_column_jacobian']={
        'weighting':'each broad source family (48-atom bulk parent vs 114-atom derived slab) has equal total squared row weight; bulk parent perturbations remain one family',
        'column_order':[c for c,a in zip(cols,active) if a],
        'singular_values':sfb.tolist(),'relative_singular_values':[float(x/sfb[0]) for x in sfb],
        'condition_number':float(sfb[0]/sfb[-1]),'rank_machine_tolerance':int(np.linalg.matrix_rank(Jfbu)),
        'family_component_row_counts':fam_counts}

    response=read_csv(OUT/'force_response_by_atom.csv')
    stability=[]
    for param in cols:
        rows=[r for r in response if r['parameter']==param]
        delta2=curv2=0.0
        for r in rows:
            for comp in 'xyz':
                b=float(r[f'base_f{comp}']);m=float(r[f'minus1_f{comp}']);p=float(r[f'plus1_f{comp}'])
                delta2+=(p-m)**2;curv2+=(p+m-2*b)**2
        ratio=math.sqrt(curv2)/(math.sqrt(delta2) or 1e-300)
        stability.append({'parameter':param,'curvature_to_central_span':ratio,
                          'status':'ZERO_RESPONSE' if delta2==0 else ('LOW_CURVATURE' if ratio<0.02 else 'NONLINEARITY_REVIEW')})
    analysis['finite_difference_stability']={x['parameter']:{'curvature_to_central_span':x['curvature_to_central_span'],'status':x['status']} for x in stability}
    (OUT/'jacobian_analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')
    with (OUT/'finite_difference_stability.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['parameter','curvature_to_central_span','status'],lineterminator='\n');w.writeheader();w.writerows(stability)

    # Periodic minimum Si-C distances check whether the Si-C radial pair enters
    # LAMMPS's Tersoff taper (R+D) in any available frozen geometry.
    label_manifest=json.loads((E4/'dft_label_manifest.json').read_text())
    grouped={}
    for r in read_csv(E4/'dft_force_labels.csv'):
        grouped.setdefault(r['case'],[]).append(r)
    si_c=[]
    for meta in label_manifest:
        cell=np.asarray(meta['cell_A'],float)
        si=[np.array([float(a[k]) for k in ('x_A','y_A','z_A')]) for a in grouped[meta['case']] if a['element']=='Si']
        cc=[np.array([float(a[k]) for k in ('x_A','y_A','z_A')]) for a in grouped[meta['case']] if a['element']=='C']
        best=math.inf
        for x in si:
            for y in cc:
                d=y-x
                for i in (-1,0,1):
                    for j in (-1,0,1):
                        for k in (-1,0,1): best=min(best,float(np.linalg.norm(d+np.array([i,j,k])@cell)))
        si_c.append({'case':meta['case'],'minimum_Si_C_distance_A':best,'pair_cutoff_R_A':2.960379,'cutoff_width_D_A':0.251126,'outer_taper_R_plus_D_A':3.211505,'inside_taper':best<=3.211505})
    with (OUT/'pair_activation_analysis.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(si_c[0]),lineterminator='\n');w.writeheader();w.writerows(si_c)
    analysis['Si_C_pair_activation']={'outer_taper_A':3.211505,'min_distance_across_cases_A':min(x['minimum_Si_C_distance_A'] for x in si_c),'all_cases_outside_taper':all(not x['inside_taper'] for x in si_c)}
    (OUT/'jacobian_analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')

    corrs=read_csv(OUT/'parameter_correlation.csv')
    pair=[r for r in corrs if r['parameter_i']<r['parameter_j'] and r['response_cosine'] not in ('', 'None')]
    pair.sort(key=lambda r:abs(float(r['response_cosine'])),reverse=True)
    elem=read_csv(OUT/'force_response_by_element.csv')
    agg={}
    for r in elem:
        k=(r['parameter'],r['element']);a=agg.setdefault(k,[0.0,0,0.0])
        n=int(r['atom_count']); rms=float(r['response_rms_vector_eV_A_per_ln_p'])
        a[0]+=n*rms*rms;a[1]+=n;a[2]+=float(r['dot_force_residual_J_half_dSSE'])
    corr_table='\n'.join(f"| {r['parameter_i']} | {r['parameter_j']} | {float(r['response_cosine']):.6f} |" for r in pair[:8])
    element_table='\n'.join(f'| {p} | {e} | {math.sqrt(v[0]/v[1]):.5g} | {v[2]:.5g} |' for (p,e),v in agg.items())
    distance_table='\n'.join(f"| {x['case']} | {x['minimum_Si_C_distance_A']:.4f} | {x['outer_taper_R_plus_D_A']:.4f} | {x['inside_taper']} |" for x in si_c)
    s=analysis['raw_jacobian']['singular_values']
    active=analysis['unit_column_jacobian_active_columns_only']
    active_cols=active['column_order']
    candidate=None;best_cond=math.inf
    for i,a in enumerate(active_cols):
        for b in active_cols[i+1:]:
            ia,ib=cols.index(a),cols.index(b)
            M=J[:,[ia,ib]]
            nrm=np.linalg.norm(M,axis=0)
            sv=np.linalg.svd(M/nrm,compute_uv=False)
            cond=float(sv[0]/sv[-1])
            if cond<best_cond:best_cond=cond;candidate=(a,b)
    E4_commit=json.loads((E4/'E4_FINAL_STATUS.json').read_text())
    mapping_sha=next((line.split()[0] for line in (E3/'SHA256SUMS.txt').read_text().splitlines() if 'TERSOFF_PARAMETER_MAPPING_AND_BLOCKS.md' in line),'not_found')

    (OUT/'PARAMETER_BLOCK_MAPPING.md').write_text(f'''# Stage E5 tied-block mapping audit

The E3 map identifies `A`, `B`, `lambda1`, and `lambda2` as equal-valued radial fields shared by each unordered first/second species pair across the third-element triplet slots. E4's mapping audit also checked the matching LAMMPS 22 Jul 2025 Update 4 `pair_tersoff.cpp` implementation (SHA-256 `5edfb3b611e3cc6d07f7cc894eb7354475b607f5e0db87d510a0bd62ceb08d93`): radial pair terms read `[i,j,j]`, while the angular/bond-order paths use triplet entries. Accordingly, this test scales all six ordered triplet rows in one unordered-pair field block together. No individual triplet row was perturbed.

Original TSCBOP SHA-256: `{analysis['source_tscbop_sha256']}`. E4 source map SHA-256: `{mapping_sha}`. The original file was copied into `potentials/TSCBOP.source.tersoff`; the source path was only read. E3 reports that paper parameters were externally compared, while E-local access to the paper was not verified; this audit does not assert independent paper verification or establish fit bounds.

| Block | Field | Pair | Base value | Tied table rows | −1% value | +1% value |
|---|---|---|---:|---:|---:|---:|
'''+''.join(f"| {p} | {field} | {a}-{b} | {float(next(v['base_value'] for v in read_csv(OUT/'parameter_variants.csv') if v['parameter']==p and v['sign']=='minus1')):.12g} | 6 | {float(next(v['perturbed_value'] for v in read_csv(OUT/'parameter_variants.csv') if v['parameter']==p and v['sign']=='minus1')):.12g} | {float(next(v['perturbed_value'] for v in read_csv(OUT/'parameter_variants.csv') if v['parameter']==p and v['sign']=='plus1')):.12g} |\n" for p,a,b,field in [('A_Ti-C','Ti','C','A'),('B_Ti-C','Ti','C','B'),('lambda1_Ti-C','Ti','C','lambda1'),('lambda2_Ti-C','Ti','C','lambda2'),('A_Ti-Si','Ti','Si','A'),('B_Ti-Si','Ti','Si','B'),('A_Si-C','Si','C','A'),('B_Si-C','Si','C','B')])+f'''

Each variant SHA and changed-row count is in `parameter_variants.csv`; `parameter_blocks.csv` records the six line/triplet members and values. The 16 variants each start from the same source file and change only the selected field in its six tied rows. All other Tersoff fields remain unchanged.
''')

    (OUT/'RUN0_REPRODUCTION_REPORT.md').write_text(f'''# Stage E5 frozen run-0 audit

The E4 package was checked against its 184-entry SHA-256 manifest. It contains six converged force-label geometries totaling 354 atoms: five 48-atom configurations from one bulk parent and one 114-atom interface-derived TSC-only slab. Actual force arrays, structure IDs, elements, coordinates and source geometry hashes were preserved. No force labels were inferred from RMSE.

LAMMPS {analysis['lammps_version']} was run with `units metal` (forces eV/Å), `pair_style tersoff`, the same `Ti Si C` map, `p p p`, neighbor settings, 1 OpenMP thread, and `run 0 post no`. The 102 calls are six common baseline checks plus 8 blocks × 2 signs × 6 structures. There was no integrator or minimizer. The frozen `atoms.data` SHA is common across baseline and all variants for each structure; all atom IDs and periodic coordinates were validated. Maximum wrapped-coordinate discrepancy was `1.07e-14 Å`.

The new baseline forces match the E4 baseline exactly for every atom/component in all six structures (`baseline_reproduction.json`). The original potential SHA remains `{analysis['source_tscbop_sha256']}`. LAMMPS-only startup diagnostics: the default singleton MPICH/UCX transport could not create a sandbox-denied local socket and exited before any force output. The successful runs used `UCX_TLS=self`, which is sufficient for the single-process executable; the scientific input/settings were otherwise unchanged. The failed startup log is retained under `run0/`.

Every successful run has its own `run.in`, frozen `atoms.data`, force dump/CSV, screen/log/launcher output and `run_metadata.json` under `run0_ucx_self/`. All records report step 0; no DFT, fitting, minimization, MD or TTM was run.
''')

    (OUT/'SVD_ANALYSIS.md').write_text(f'''# Stage E5 force-response Jacobian and SVD

The matrix has **{J.shape[0]} force-component rows × 8 tied-block columns**, assembled from 354 atom vectors in fixed case/atom/component order. Each column is the central finite difference `J=[F(1.01p)−F(0.99p)]/0.02`, i.e. an approximate derivative with respect to `ln(p)`, not an optimized parameter direction.

## Rank and conditioning

- Full 8-column Jacobian numerical rank: **{analysis['raw_jacobian']['rank_machine_tolerance']}/8**. The smallest two singular values are exactly zero because both Si–C columns have no force response in this dataset; its condition number is infinite/undefined.
- The remaining six columns have rank **6/6** at machine tolerance. After unit-normalizing active columns, condition number is **{active['condition_number']:.3f}**; relative singular values: `{', '.join(f'{x:.4g}' for x in active['relative_singular_values'])}`. The effective rank is 5 at relative threshold `1e-3` and 4 at `1e-2`, showing poor separation despite full algebraic rank.
- A family-balanced sensitivity check (bulk-parent and derived-slab groups given equal total row weight) has active-column condition number **{analysis['family_balanced_unit_column_jacobian']['condition_number']:.3f}**. The conclusion remains ill-conditioned and two Si–C columns stay inactive.

## Highest response-column correlations

| Parameter i | Parameter j | Cosine similarity |
|---|---|---:|
{corr_table}

Most important coupled directions include `B_Ti-C`/`lambda1_Ti-C` (+0.999881), `A_Ti-C`/`lambda1_Ti-C` (−0.998078), and `A_Ti-Si`/`B_Ti-Si` (−0.998186). These finite-difference force columns cannot robustly identify all active radial fields together from the present geometries.

## Element-wise response magnitude and local residual direction

RMS response is vector RMS over atoms of each element, in eV/Å per unit `ln(p)`. The final column is the summed baseline force-residual dot response over six structures; its sign indicates only a local direction for squared force error, not a fitted result.

| Parameter | Element | RMS response | Σ residual·J |
|---|---|---:|---:|
{element_table}

The five bulk cases share one parent structure, and the 114-atom slab is a derived diagnostic structure. Their agreement is a consistency check, not independent validation. In particular, the Ti–Si parameter directions for Si residuals conflict between the bulk-parent and slab groups (`structure_family_consistency.csv`). `lambda2_Ti-C` also has opposite residual-direction signs for Ti and C, so the global scalar error direction hides an element tradeoff.

Files: `PARAMETER_RESPONSE_MATRIX.csv`, `force_response_by_atom.csv`, `force_response_by_element.csv`, `force_response_by_structure.csv`, `parameter_correlation.csv`, `structure_family_consistency.csv`, `finite_difference_stability.csv`, and `jacobian_analysis.json`.
''')

    inactive=[c for c,n in zip(cols,np.linalg.norm(J,axis=0)) if n<=max(np.linalg.norm(J,axis=0))*1e-12]
    chosen=', '.join(candidate) if candidate else 'none'
    (OUT/'IDENTIFIABLE_PARAMETER_SET.md').write_text(f'''# Stage E5 identifiable parameter response set

**Observed force-active columns:** `{', '.join(active_cols)}`. They form a full-rank six-column matrix algebraically, but the unit-column condition number is {active['condition_number']:.1f}; singular directions are weak and several radial columns are nearly collinear. In particular, do not simultaneously free the four Ti–C radial fields or both Ti–Si amplitude fields from this dataset.

**Unexcited columns:** `{', '.join(inactive)}` have exactly zero response in every tested atom/component. The minimum periodic Si–C distances are all outside the Tersoff Si–C taper, whose outer limit is `R+D=3.211505 Å`:

| Structure | Minimum Si–C distance (Å) | Outer taper (Å) | Within taper |
|---|---:|---:|---|
{distance_table}

Therefore the present force labels cannot constrain `A_Si-C` or `B_Si-C`; those pair blocks need geometries that activate Si–C interactions before any fit involving them. This zero response is a structural coverage gap, not evidence that the parameters are physically irrelevant.

Among two-column subsets of the six active responses, the best-conditioned pair in this local matrix is `{chosen}` (unit-column condition number {best_cond:.3f}). Treat this only as a candidate pair for a future, explicitly approved exploratory test: data-family independence, correct equilibrium/elastic response and validation roles are still missing. It does not establish unique parameter values or generalization.
''')

    (OUT/'TRIAL_FIT_RECOMMENDATION.md').write_text(f'''# Stage E5 trial-fit gate

**Recommendation: START_TRIAL_FIT=NO.** This stage completed the authorized static force-response analysis only. The Jacobian shows two unexcited parameters, high coupling among active fields, and no independent DEV/TEST family. The current 48-atom perturbations all belong to one parent; the 114-atom TSC-only slab is derived and is not an independent clean-surface validation set. No roles were changed.

The smallest numerically separated two-column candidate is `{chosen}` (condition {best_cond:.3f} after unit-column scaling), but the response matrix is not a substitute for a balanced DFT objective or held-out validation. Before any separately authorized fit, obtain Si–C-active geometries, independent volume/shear/Si-displacement and clean-surface data as E3 specified, then freeze TRAIN/DEV/TEST by source family. If a fit is later approved, hold all other Tersoff blocks fixed and regularize toward the original parameterization; the ±1% probe itself is not a bound.

No TSCBOP, EAM, Morse, DFT labels, or frozen structures were modified. No fitting/optimization, DFT, MD, TTM, or MACE run was started. Candidate parameters remain diagnostic only and are not production-ready.
''')

    (OUT/'E5_FINAL_REPORT.md').write_text(f'''# Stage E5 final report — Ti₃SiC₂ TSCBOP force-response Jacobian

## Findings

- Eight tied scalar blocks were tested separately at −1%/+1%; 102 single-process LAMMPS `run 0` evaluations include common baselines. Each block changes only the six rows sharing its unordered pair/field.
- E4 DFT-force labels were hash-verified and contain 354 atoms; all six new baselines reproduce E4 forces exactly. The source TSCBOP SHA remains `{analysis['source_tscbop_sha256']}`.
- Full Jacobian shape: `{J.shape[0]}×{J.shape[1]}`; rank **{analysis['raw_jacobian']['rank_machine_tolerance']}/8**. `A_Si-C` and `B_Si-C` are zero-response columns because every available Si–C distance exceeds their outer cutoff. Active six-column rank is 6, but normalized condition is {active['condition_number']:.1f}; effective rank falls to 5 at 1e-3 and 4 at 1e-2.
- Strong coupling prevents robust simultaneous use of the four Ti–C radial blocks and Ti–Si A/B. The best-conditioned active two-column slice is `{chosen}` (condition {best_cond:.3f}); this is a local response-screening result only.
- The five 48-atom cases are one parent family; the one 114-atom slab is derived. Thus there is no independent validation family and no claim of general parameter identifiability or predictive performance.

## Required decision fields

| Field | Result | Basis |
|---|---|---|
| `PARAMETER_BLOCK_MAPPING_PASS` | YES | E3 shared-pair map, original file, exact tied-row edits, and LAMMPS implementation mapping |
| `RAW_FORCE_DATA_AVAILABLE` | YES | Six converged E4 geometries, actual hashed force labels |
| `RUN0_REPRODUCED` | YES | All six E5 source baselines match E4 per component exactly |
| `FORCE_RESPONSE_COMPUTED` | YES | Eight columns, 354 atom vectors, 1062 force components |
| `JACOBIAN_IDENTIFIABILITY` | PARTIAL / ILL-CONDITIONED | Six active algebraic directions, two inactive columns, high correlations and no independent validation family |
| `IDENTIFIABLE_PARAMETER_SET` | `{', '.join(active_cols)} (local force-active only); Si–C A/B not constrained` | Rank and distance/taper audit |
| `TRIAL_FIT_RECOMMENDATION` | NO | Leave all fitting on hold pending data coverage and independent family roles |
| `PRODUCTION_READY` | NO | No fit and no independent validation |

The output package includes raw per-atom forces for every baseline/variant, parameter response matrix, per-element/per-structure summaries, correlations, SVD analysis, finite-difference stability, pair activation checks, and SHA-256 manifest. No fit, DFT, MD, TTM, or MACE training was run.
''')

    print('E5 reports written; active columns:',', '.join(active_cols),'inactive:',', '.join(inactive),'best two:',chosen,'condition:',best_cond)

if __name__=='__main__': main()
