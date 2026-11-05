# Open-top checkpoint restart smoke

The open-top 0.160 ps checkpoint was resumed using matching pressure, source-ramp, and boundary settings. The solver restored 4130.42566868 nm³ Ag and 1784.17049967 nm³ skeleton and reached 0.161 ps after 18 steps (exit code 0). A compressed 0.161 ps checkpoint and sidecar are retained here. To pass the compressed dump back to the current executable, decompress it to the matching `.dump` path and keep its `.state` and `.complete` sidecars alongside it.
