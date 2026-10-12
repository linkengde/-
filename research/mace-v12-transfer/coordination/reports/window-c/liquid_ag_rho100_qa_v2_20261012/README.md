# D Ag64 rho100 QA V2 evidence package

This package contains the read-only QA V2 implementation, deterministic synthetic tests, environment fingerprint, and regression evidence for C review. It does not contain the Ag EAM potential or a runnable MD continuation runner. The original step-1000 checkpoint and historical trajectory are referenced by path and SHA only; neither was modified.

The package passes its engineering tests, but the archived local-Q6 regression misses the C lock's documented target. The lock remains unsigned and the config keeps `execution_authorized=false`. Do not start or resume MD from this package. See `D_QA_V2_ADAPTATION_REPORT.md` and `Q6_REGRESSION_REPORT.json`.
