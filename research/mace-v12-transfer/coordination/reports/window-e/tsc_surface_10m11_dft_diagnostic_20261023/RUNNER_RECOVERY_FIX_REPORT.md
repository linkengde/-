# SCF result-reader recovery

The GPAW log for the resumed run records `Converged in 24 steps` at global SCF
iteration 64 and lists all three convergence indicators as passed. The process
then exited with code 42 because the original collector required both that log
marker and an ASE `converged` hook flag. On this GPAW 26.7 restart path, the hook
flag was false despite the authoritative convergence log. No final label file
was written; no intermediate values are promoted to labels.

A separate recovery runner was derived from the original runner. Its only
behavioral change is to require the GPAW log convergence marker without making
the hook flag a second gate. It records both the log marker and hook value. The
original runner, input, scientific config, and original execution lock remain
unchanged. The recovery lock binds the derived runner SHA and retains the same
configuration, geometry, input, PAW, and scientific-parameter hashes.

Recovery source: verified `checkpoint_iter_0060.gpw` (SHA-256
`5cd9f365ac2b80e888e3f357cb93a12798f59c263cd6b4b823017c16a60d6c2a`), read-only
4-rank reload PASS with no SCF. The recovery recomputes the final SCF steps and
collects labels only after GPAW again reports convergence.
