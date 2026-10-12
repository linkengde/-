# Previous R1 attempt retained as audit evidence

A run had begun under the immediately preceding explicit run instruction. A newer user instruction then required submitting an execution authorization request before any DFT execution. The MPI ranks were terminated at iteration 2.

- State: `ABORTED_PENDING_C_AUTH`
- SCF converged: no
- Completed iterations: 2
- Valid energy/force/stress result: no
- Checkpoint: none
- The partial GPAW and launcher logs are audit evidence only and must not be used as a DFT label, ranking input, or training data.
- No R1/R2 comparison or ranking was updated from that partial attempt.

The next run must use a fresh independent output directory and must pass the formal C authorization gate before GPAW starts.
