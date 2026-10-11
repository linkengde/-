# Reference impact on Ag–Si DEV v2

Abdul Wasey et al. (2013), DOI `10.1002/pssb.201248542`, is now verified from the user-provided primary PDF. It supports using coincidence-site epitaxy as a design method, considering a 4×Ag(111)/3×Si(111) match, examining the paper's illustrated registry (every fifth Ag atom above the fourth Si atom), and using its reported local Ag–Si bond lengths (about 2.53–2.72 Å) as qualitative scale information. The paper also documents a slab-based VASP/PAW relaxation workflow.

Its impact on the proposed Ag/Ti₃SiC₂ Si-termination DEV is limited to **design-method guidance**. Elemental Si(111) is not the Ti₃SiC₂(0001) surface: the actual MAX-phase in-plane lattice, termination motif, strain, registry, and interface spacing must be constructed and verified independently. In particular, the paper reports local Ag–Si bond lengths, not a plane-to-plane (d_0); its single illustrated registry is not proof of a global registry minimum.

The article does not provide a complete machine-readable geometry or per-atom forces for the target Ag/Ti₃SiC₂ interface. Its values therefore cannot be used as Ag/Ti₃SiC₂ training labels, and the paper's VASP/LDA/PBE settings are not automatically the project's GPAW/PBE input settings.

This record does **not**:

- copy an Ag/Si geometry into the Ag/Ti₃SiC₂ data set;
- create or relabel TRAIN, DEV, or TEST data;
- establish source-family independence or clear the previous target-structure blockers;
- authorize structure generation, DFT, or MACE training.

The target Si-terminated Ti₃SiC₂ slab, layer/atom map, in-plane matching, explicit registry coordinates, target (d_0), vacuum/PBC, and hashed input remain unresolved. Existing data roles and numerical-control labels remain unchanged.
