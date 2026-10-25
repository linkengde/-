#!/usr/bin/env python3
"""Build the CIF-derived 3x3x2 orthorhombic Ti3SiC2 TSCBOP probe cell."""
from pathlib import Path
from ase.build import make_supercell
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[4]
CIF = ROOT / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
OUTPUT = Path(__file__).resolve().parent / "ti3sic2_bulk_3x3x2.data"

atoms = read(CIF)
atoms = make_supercell(atoms, [[1, 0, 0], [1, 2, 0], [0, 0, 1]])
atoms = atoms.repeat((3, 3, 2))
atoms.pbc = True
write(OUTPUT, atoms, format="lammps-data", atom_style="atomic", specorder=["Ti", "Si", "C"], masses=True)
print(f"Wrote {len(atoms)} atoms, formula {atoms.get_chemical_formula()}, cell lengths {atoms.cell.lengths()} A")
