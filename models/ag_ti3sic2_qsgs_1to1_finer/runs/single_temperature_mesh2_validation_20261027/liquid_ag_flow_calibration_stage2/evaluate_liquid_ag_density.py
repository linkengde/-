#!/usr/bin/env python3
"""Evaluate the source-bounded Abdullaev et al. (2020) liquid-Ag density fit.

This utility evaluates a published experimental correlation; it does not
produce raw experimental measurements or extrapolate outside the source range.
"""

import argparse
import json


T_MIN_K = 1234.93
T_MAX_K = 1631.0
RHO0_KG_M3 = 9327.1
SLOPE_KG_M3_K = 0.9956
SOURCE = (
    "R. N. Abdullaev, R. A. Khairulin, S. V. Stankus, "
    "J. Phys.: Conf. Ser. 1677 (2020) 012161, Eq. (2), "
    "DOI: 10.1088/1742-6596/1677/1/012161"
)


def density_kg_m3(temperature_k: float) -> float:
    if not T_MIN_K <= temperature_k <= T_MAX_K:
        raise ValueError(
            f"T={temperature_k:g} K is outside the source range "
            f"[{T_MIN_K:g}, {T_MAX_K:g}] K; extrapolation is disabled"
        )
    return RHO0_KG_M3 - SLOPE_KG_M3_K * (temperature_k - T_MIN_K)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--temperature-k",
        type=float,
        action="append",
        dest="temperatures",
        help="Temperature in K; repeat to evaluate multiple points",
    )
    args = parser.parse_args()
    temperatures = args.temperatures or [T_MIN_K, 1250.0, 1300.0, 1600.0, T_MAX_K]
    rows = [
        {"temperature_K": t, "density_kg_m3": density_kg_m3(t)}
        for t in temperatures
    ]
    print(
        json.dumps(
            {
                "status": "SOURCE_BOUNDED_CORRELATION_EVALUATION",
                "equation": "rho=9327.1-0.9956*(T-1234.93)",
                "temperature_range_K": [T_MIN_K, T_MAX_K],
                "source": SOURCE,
                "uncertainty_note": (
                    "Paper reports density errors 0.25% at Tm, 0.26% at "
                    "1300 K and 0.34% at 1600 K; coverage convention is unspecified. "
                    "No pointwise 1250 K uncertainty is inferred."
                ),
                "values": rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
