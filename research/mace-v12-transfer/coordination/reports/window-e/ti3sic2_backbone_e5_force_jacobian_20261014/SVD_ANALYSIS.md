# Stage E5 force-response Jacobian and SVD

The matrix has **1062 force-component rows × 8 tied-block columns**, assembled from 354 atom vectors in fixed case/atom/component order. Each column is the central finite difference `J=[F(1.01p)−F(0.99p)]/0.02`, i.e. an approximate derivative with respect to `ln(p)`, not an optimized parameter direction.

## Rank and conditioning

- Full 8-column Jacobian numerical rank: **6/8**. The smallest two singular values are exactly zero because both Si–C columns have no force response in this dataset; its condition number is infinite/undefined.
- The remaining six columns have rank **6/6** at machine tolerance. After unit-normalizing active columns, condition number is **1104.743**; relative singular values: `1, 0.2702, 0.05047, 0.01742, 0.003777, 0.0009052`. The effective rank is 5 at relative threshold `1e-3` and 4 at `1e-2`, showing poor separation despite full algebraic rank.
- A family-balanced sensitivity check (bulk-parent and derived-slab groups given equal total row weight) has active-column condition number **903.158**. The conclusion remains ill-conditioned and two Si–C columns stay inactive.

## Highest response-column correlations

| Parameter i | Parameter j | Cosine similarity |
|---|---|---:|
| B_Ti-C | lambda1_Ti-C | 0.999881 |
| A_Ti-Si | B_Ti-Si | -0.998186 |
| A_Ti-C | lambda1_Ti-C | -0.998078 |
| A_Ti-C | B_Ti-C | -0.997647 |
| B_Ti-C | lambda2_Ti-C | -0.994693 |
| lambda1_Ti-C | lambda2_Ti-C | -0.994037 |
| A_Ti-C | lambda2_Ti-C | 0.985477 |
| A_Ti-Si | lambda2_Ti-C | -0.857104 |

Most important coupled directions include `B_Ti-C`/`lambda1_Ti-C` (+0.999881), `A_Ti-C`/`lambda1_Ti-C` (−0.998078), and `A_Ti-Si`/`B_Ti-Si` (−0.998186). These finite-difference force columns cannot robustly identify all active radial fields together from the present geometries.

## Element-wise response magnitude and local residual direction

RMS response is vector RMS over atoms of each element, in eV/Å per unit `ln(p)`. The final column is the summed baseline force-residual dot response over six structures; its sign indicates only a local direction for squared force error, not a fitted result.

| Parameter | Element | RMS response | Σ residual·J |
|---|---|---:|---:|
| A_Ti-C | Ti | 26.142 | 10411 |
| A_Ti-C | Si | 0 | 0 |
| A_Ti-C | C | 2.8933 | 751.44 |
| B_Ti-C | Ti | 25.333 | -10076 |
| B_Ti-C | Si | 0.44934 | -30.331 |
| B_Ti-C | C | 0.8732 | -190.11 |
| lambda1_Ti-C | Ti | 79.603 | -31701 |
| lambda1_Ti-C | Si | 0 | 0 |
| lambda1_Ti-C | C | 2.989 | -741.97 |
| lambda2_Ti-C | Ti | 60.218 | 23945 |
| lambda2_Ti-C | Si | 1.5385 | 103.95 |
| lambda2_Ti-C | C | 5.5442 | -1439.1 |
| A_Ti-Si | Ti | 6.5321 | -2519.4 |
| A_Ti-Si | Si | 4.6371 | -280.28 |
| A_Ti-Si | C | 0 | 0 |
| B_Ti-Si | Ti | 7.7974 | 3033.4 |
| B_Ti-Si | Si | 5.984 | 369.25 |
| B_Ti-Si | C | 0.041666 | -9.0946 |
| A_Si-C | Ti | 0 | 0 |
| A_Si-C | Si | 0 | 0 |
| A_Si-C | C | 0 | 0 |
| B_Si-C | Ti | 0 | 0 |
| B_Si-C | Si | 0 | 0 |
| B_Si-C | C | 0 | 0 |

The five bulk cases share one parent structure, and the 114-atom slab is a derived diagnostic structure. Their agreement is a consistency check, not independent validation. In particular, the Ti–Si parameter directions for Si residuals conflict between the bulk-parent and slab groups (`structure_family_consistency.csv`). `lambda2_Ti-C` also has opposite residual-direction signs for Ti and C, so the global scalar error direction hides an element tradeoff.

The curvature-to-central-span ratio for each ±1% probe is at most 0.0136 (`finite_difference_stability.csv`), so the tested central differences are locally stable at this step size. This does not establish behavior outside the ±1% neighborhood.

Files: `PARAMETER_RESPONSE_MATRIX.csv`, `force_response_by_atom.csv`, `force_response_by_element.csv`, `force_response_by_structure.csv`, `parameter_correlation.csv`, `structure_family_consistency.csv`, `finite_difference_stability.csv`, and `jacobian_analysis.json`.
