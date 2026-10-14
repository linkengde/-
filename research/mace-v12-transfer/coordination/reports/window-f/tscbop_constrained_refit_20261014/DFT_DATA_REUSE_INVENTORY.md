# DFT data reuse inventory

Audit performed on branch `cloud-f-tscbop-refit`, based on E commit `7701ecfdd4e565c4d4be16b889d091c44aa8b648`. Only the explicitly authorized pure-phase Ti₃SiC₂ records and the published E reports were inspected. No sealed test data were opened.

## Available 48-atom bulk records

All seven result directories pass their local `SHA256SUMS.txt`. The three baseline meshes share input SHA-256 `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`, exact ordered species/positions/cell, and canonical geometry SHA-256 `269c1be824ccd905bb7ab05bf8e955f3f9d2d9c77839d2ff95fef2d074f5d05b`. They are one geometry with numerical mesh repeats, not three independent structures. The result archives identify the role as numerical pure-phase controls, not training or independent validation.

| Record | k mesh | Result extxyz SHA-256 | DFT E (eV/cell) | Free E (eV/cell) | SCF | fmax (eV/Å) | Geometry/role |
|---|---:|---|---:|---:|---|---:|---|
| baseline k4 | 4×4×2 | `20f15465deb0b6f1661385e632ff4a51686fecc4654bba84bc11f08d898e691c` | −389.9010449657 | −390.1250638853 | converged, 46 iterations | 0.09424 | shared parent geometry |
| baseline k6 | 6×6×2 | `fe7ec9ce8b17884aab37f45d2c1c494878c7536571ca41a41bbf8c584e309c8c` | −389.8177237814 | −390.0641824804 | converged, 41 iterations | 0.08418 | same geometry; mesh repeat |
| baseline k10 | 10×10×4 | `75b083e225b8b20bbdde181c8e1847bae4b0052a854c6e9ee48885c940166415` | −389.8092428601 | −390.0617903646 | converged, 34 iterations | 0.08786 | same geometry; mesh repeat |

The records use GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, PBE, 500 eV plane waves, 0.10 eV smearing, 4 MPI ranks, and native extrapolated energy with free energy recorded separately. The baseline source geometry is CIF-derived; the archive explicitly limits it to a numerical control. It is not a relaxed equilibrium cell or an independent validation structure. The extxyz `Properties` record stores species, positions, forces, and local IDs; no stress property is present.

## Existing Ti/C displacement responses

Each displacement result is converged at 10×10×4 with the same DFT method and same parent cell. The two Ti cases move only local ID 3 in z by ±0.020 Å; the two C cases move only local ID 9 in z by ±0.020 Å. Their parent proposal SHA-256 values and result hashes are:

| Record | Proposal SHA-256 | Result extxyz SHA-256 | E (eV/cell) | Free E (eV/cell) | fmax (eV/Å) | Role/limit |
|---|---|---|---:|---:|---:|---|
| Ti ID 3, +0.020 Å z | `5829ebcbb9ae331f6cf1dd7dc5abca0787c07667a55778e3c1eea241938b227b` | `fa16c172aa4ae926b38347c230e60f8bb93dd79a2e2d66ec2481de79579d1265` | −389.8045908244 | −390.0574330057 | 0.39058 | same parent family |
| Ti ID 3, −0.020 Å z | `5529a040f0f35957f826b78754910970741bff741474e321dc1646a5c0c12bbf` | `95f3803ec53eaa5a89c6b0117446782e66f4060dbcb01588a7bb5488b24220ce` | −389.8069019827 | −390.0591550120 | 0.30828 | same parent family |
| C ID 9, +0.020 Å z | `6308961f33ad5fae3f58930570a47044a236e062de243821c20577c78955ff72` | `3c0d7f82a748ea93e3afc1d95da44e28f79e11e3ef3959e5746c1463c7c66229` | −389.8080388055 | −390.0605252760 | 0.21642 | same parent family |
| C ID 9, −0.020 Å z | `19f1a4011620f5c7f245d0a8c893964aa8ebd1eb647f8a9b6f4547a75cac586f` | `ac021467637026f8728b9b25228b061f8485dda42577b5ce8ae5a8ab4b134850` | −389.8043617972 | −390.0569707520 | 0.39360 | same parent family |

The local IDs and ordered species are present and unique in all four outputs. These are useful force-response points, but all four derive from the same parent geometry, sample only one crystallographic site per species and the z direction, and carry no stress labels. They cannot serve as independent validation.

## Surface and other records

The E report describes a frozen 114-atom `Ti57Si19C38` TSC-only slab derived by removing Ag from a relaxed interface. It reports hashes for a 135-atom parent (`eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f`), the pre-relaxation structure (`1a3e341e1132bba9b349ab0fe6e3a3be11d1e1661671b96f6810a54f2e1f54ae`), the 114-atom force-error table (`8028ff81383da0b151be2c2d70e5b595dc8b7f7b280cdfcf75903e9710141dd4`), and coordination table (`aac97575856d2a48d67bcc99fb5c530fefc38951f6e99845b757c2e25c08d5e4`). The source slab, raw DFT labels, convergence summary, and stress are not present in this checkout snapshot. It is therefore report-only evidence here, not a directly re-audited fitting label. The E report states its provenance as an interface-derived frozen geometry, not an independently prepared clean surface.

The pure-phase directory also contains volume and shear proposal files, but no converged DFT outputs for those proposals. They are not labels. No additional source-audited, converged pure-phase structures were verified in this audit.

## Completeness for fitting

The archived bulk force labels are enough to confirm that a large force mismatch was reported on the shared parent geometry. The available fitting data do not include stress, a volume/shear response curve, independently sourced configurations, or locally auditable raw surface labels. Thus data provenance for the inspected 48-atom records passes, while a complete fitting/validation data audit does not. No fitting was started.
