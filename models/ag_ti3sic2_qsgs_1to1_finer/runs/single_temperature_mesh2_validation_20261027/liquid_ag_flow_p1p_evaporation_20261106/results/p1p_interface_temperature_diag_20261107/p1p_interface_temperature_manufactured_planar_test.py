#!/usr/bin/env python3
"""Manufactured planar VOF-interface temperature reconstruction check.

A cell-centered temperature is the volume average of a continuous, piecewise
linear liquid/gas temperature field. The PLIC interface lies at x=q in the
cut cell [0,1]. Compare the production centered-gradient extrapolation with a
liquid-side least-squares extrapolation from three pure-liquid cells.
"""
import csv
import math
from pathlib import Path

OUT = Path("/workspace/ag70-p1p-evaporation/models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/liquid_ag_flow_p1p_evaporation_20261106/results/p1p_interface_temperature_diag_20261107/p1p_interface_temperature_manufactured_planar.csv")

def cell_average(i, q, t_surface, grad_liquid, grad_gas):
    lo, hi = float(i), float(i + 1)
    integral = 0.0
    a, b = lo, min(hi, q)
    if b > a:
        integral += t_surface*(b-a) + 0.5*grad_liquid*((b-q)**2-(a-q)**2)
    a, b = max(lo, q), hi
    if b > a:
        integral += t_surface*(b-a) + 0.5*grad_gas*((b-q)**2-(a-q)**2)
    return integral

def fit_intercept(xs, ys):
    xm, ym = sum(xs)/len(xs), sum(ys)/len(ys)
    slope = sum((x-xm)*(y-ym) for x,y in zip(xs,ys))/sum((x-xm)**2 for x in xs)
    return ym - slope*xm

rows=[]
for q in (0.1,0.3,0.5,0.7,0.9):
    for grad_liquid in (10.0,25.0,50.0):
        for ratio in (-0.5,0.0,0.5,1.0):
            t_surface=1400.0
            grad_gas=ratio*grad_liquid
            tc=cell_average(0,q,t_surface,grad_liquid,grad_gas)
            grad_centered=(cell_average(1,q,t_surface,grad_liquid,grad_gas)-
                           cell_average(-1,q,t_surface,grad_liquid,grad_gas))/2.0
            t_centered=tc+grad_centered*(q-0.5)
            xs=[i+0.5-q for i in (-1,-2,-3)]
            ys=[t_surface+grad_liquid*x for x in xs]
            t_one_sided=fit_intercept(xs,ys)
            analytical_error=0.125*(grad_gas-grad_liquid)
            assert abs((t_centered-t_surface)-analytical_error)<1e-9
            assert abs(t_one_sided-t_surface)<1e-9
            rows.append({
                "interface_fraction_q":q,
                "T_interface_exact_K":t_surface,
                "liquid_gradient_K_per_cell":grad_liquid,
                "gas_gradient_over_liquid":ratio,
                "production_centered_T_K":t_centered,
                "centered_error_K":t_centered-t_surface,
                "centered_error_analytic_K":analytical_error,
                "one_sided_liquid_fit_T_K":t_one_sided,
                "one_sided_fit_error_K":t_one_sided-t_surface,
            })
with OUT.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys())
    w.writeheader();w.writerows(rows)
max_center=max(abs(r["centered_error_K"]) for r in rows)
max_one=max(abs(r["one_sided_fit_error_K"]) for r in rows)
print(f"cases={len(rows)} max_centered_error_K={max_center:.12g} max_one_sided_error_K={max_one:.12g}")
print(f"csv={OUT}")
