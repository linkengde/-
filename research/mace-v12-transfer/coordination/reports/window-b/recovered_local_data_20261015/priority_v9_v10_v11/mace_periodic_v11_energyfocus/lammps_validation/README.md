# v11 LAMMPS gate

This folder contains the prepared two-rank `run 0` parity and one-step unheated NVE inputs. Do not execute them until v11 has been scored against v9/v10 on the frozen test and external AgSi/AgC structures. These inputs use only the 49-atom AgTi DFT validation cell and do not invoke TTM or any production model.

If the energy and force gates pass, export the v11 model with `export_v11_mliap.py`, then run the two inputs with `/workspace/.local/bin/lammps-mace-mpi-kokkos`. Save stdout/logs in this directory and compare `run 0` energy/forces against ASE using the same serialized v11 model.
