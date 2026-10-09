# Expanded-vacuum periodic control

A6/8 fixed-sigma0.10 comparison passes all numerical budgets on this geometry. Provisional recipe6x6x1/PW500/PBE/sigma0.10 retained; native extrapolated and finite-width free energies both recorded. This is not a zero-temperature or model accuracy claim.

A performs only one new numerical job: original vacuum16.00026725A ->26.00026725A, height+10A, all ions shifted+5A, same relative geometry/elements/IDs/pair and all-periodic PBC. Input is copied byte-identically from B reviewed proposal. No original baseline rerun, dipole correction, relaxation or changed occupations. Compare exact baseline label with new label using compare_mesh_results.py; it verifies translation/cell change and unchanged method. Root and result inventories preserved; four MPI ranks, one thread each,8GiB start/3GiB runtime disk reserves. Dipole/slab-off/on are separate boundary controls, not mixed into this run.
