# v10 training data

The v10 training split preserves the v9 force and interface data, adds four GPAW PW-PBE solid reference structures (fcc Ag, hcp Ti, diamond Si, diamond C), adds the DFT energy to the existing AgTi 2.3871 A frame and adds periodic AgTi 2.55 A energy/force labels. A periodic AgTi 2.70 A frame stays in the test split. The training energy composition matrix has full rank 4/4.

The pure-phase frames establish an identifiable composition basis for MACE `E0s=estimated`; they do not constitute an experimental formation-energy calibration. The 2.70 A test is a same-motif distance holdout, not independent validation of all four cross interactions or disordered/liquid interfaces. The earlier AgTi relative-energy profile remains a useful focused check.
