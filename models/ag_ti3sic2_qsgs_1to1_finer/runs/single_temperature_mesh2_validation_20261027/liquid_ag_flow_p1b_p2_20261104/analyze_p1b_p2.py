#!/usr/bin/env python3
"""Regenerate P1b grid summary and NIST Ag viscosity uncertainty tables."""
from __future__ import annotations
import csv, json, math, re
from pathlib import Path
HERE = Path(__file__).resolve().parent
OFF = HERE / "results/p1b_official"
NIST = HERE / "results/p2_nist"

def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def ag_observations(data):
    assert data["Citation"]["sDOI"] == "10.1007/s10765-016-2104-7"
    out = []
    for row in data["PureOrMixtureData"]:
        comps = row.get("Component", [])
        if not any(c.get("RegNum", {}).get("nOrgNum") == 2 for c in comps):
            continue
        prop_meta = row["Property"][0]
        assert prop_meta["Property-MethodID"]["PropertyGroup"]["TransportProp"]["ePropName"] == "Viscosity, Pa*s"
        cm = prop_meta["CombinedUncertainty"]
        for nv in row["NumValues"]:
            T = float(next(v["nVarValue"] for v in nv["VariableValue"] if v["nVarNumber"] == 1))
            pv = next(v for v in nv["PropertyValue"] if v["nPropNumber"] == 1)
            out.append((T, float(pv["nPropValue"]), float(pv["CombinedUncertainty"]["nCombExpandUncertValue"])))
        assert cm["nCombUncertLevOfConfid"] == 95
        assert cm["eCombUncertEvalMethod"] == "Propagation of evaluated standard uncertainties"
    return sorted(out)

raw_path = NIST / "nist_original_10.1007_s10765-016-2104-7.json"
sub_path = NIST / "nist_ag_9points_extracted_not_original.json"
raw = json.loads(raw_path.read_text()); sub = json.loads(sub_path.read_text())
raw_obs = ag_observations(raw); sub_obs = ag_observations(sub)
assert len(raw_obs) == len(sub_obs) == 9
assert raw_obs == sub_obs, "user-provided extracted Ag subset differs from original NIST JSON"

def fit(T):
    return math.exp(-0.4631 + 2089.0/T) # mPa s

nist_rows=[]; fit_rows=[]
for T, mu_pa, U_pa in raw_obs:
    mu=1000*mu_pa; U=1000*U_pa; pred=fit(T)
    nist_rows.append({"temperature_K":f"{T:.0f}","viscosity_Pa_s":f"{mu_pa:.9g}","expanded_uncertainty_95_Pa_s":f"{U_pa:.9g}","viscosity_mPa_s":f"{mu:.9g}","expanded_uncertainty_95_mPa_s":f"{U:.9g}","relative_expanded_uncertainty_pct":f"{100*U/mu:.6f}"})
    fit_rows.append({"temperature_K":f"{T:.0f}","observed_mPa_s":f"{mu:.9g}","expanded_uncertainty_95_mPa_s":f"{U:.9g}","arrhenius_fit_mPa_s":f"{pred:.9g}","fit_minus_observed_mPa_s":f"{pred-mu:.9g}","relative_fit_deviation_pct":f"{100*(pred-mu)/mu:.6f}","residual_over_expanded_uncertainty":f"{(pred-mu)/U:.6f}","status":"literature correlation evaluated within NIST observation range"})
write_csv(NIST/"ag_viscosity_experimental.csv", nist_rows)
write_csv(NIST/"ag_viscosity_fit_vs_data.csv", fit_rows)

# Parse official Basilisk output. Each 'logfile' section starts at its header.
stdout=(OFF/"official_stdout.log").read_text().splitlines()
blocks=[]; current=None
for line in stdout:
    if line.startswith("t dsb "):
        if current is not None: blocks.append(current)
        current=[]
    elif current is not None:
        fields=line.split()
        try:
            vals=[float(x) for x in fields[:7]]
        except (ValueError, IndexError):
            continue
        if len(vals)==7: current.append(vals)
if current is not None: blocks.append(current)
assert len(blocks)==3, f"expected 3 resolution blocks, got {len(blocks)}"
mesh=[]
for cells, rows in zip((8,16,32),blocks):
    stderr=(OFF/"official_stderr.log").read_text().splitlines()
    m=next((re.match(rf"^{cells} ([0-9.eE+-]+) ([0-9.eE+-]+)$",s.strip()) for s in stderr if re.match(rf"^{cells} ([0-9.eE+-]+) ([0-9.eE+-]+)$",s.strip())),None)
    assert m, f"missing final speed for R/dx={cells}"
    measured_over_U=float(m.group(1)); theory_over_U=float(m.group(2))
    ratio=measured_over_U/theory_over_U
    mass_pct=[100*r[1] for r in rows]
    mesh.append({"cells_per_radius":cells,"N_domain":cells*16,"final_sample_t_over_t0":f"{rows[-1][0]:.9g}","measured_velocity_over_U":f"{measured_over_U:.9g}","theory_velocity_over_U":f"{theory_over_U:.9g}","measured_over_theory":f"{ratio:.9g}","relative_velocity_error_pct":f"{100*abs(ratio-1):.6f}","final_logged_mass_drift_pct":f"{mass_pct[-1]:.9g}","max_abs_logged_mass_drift_pct":f"{max(abs(x) for x in mass_pct):.9g}","log_samples":len(rows),"pass_speed_trend":""})
for i,row in enumerate(mesh):
    row["pass_speed_trend"] = "coarsest baseline" if i==0 else ("error decreases with refinement" if float(row["relative_velocity_error_pct"]) < float(mesh[i-1]["relative_velocity_error_pct"]) else "error did not decrease")
write_csv(OFF/"p1b_velocity_grid_convergence.csv",mesh)

# The project's automatic-curvature VOF comparison uses the same physical case
# but a separate standard-VOF/CSF implementation, not the official CLSVOF path.
AUTO = HERE / "results/p1b_autocurvature"
auto_stdout=(AUTO/"p1b_vof_fullstress_stdout.log").read_text().splitlines()
auto_stderr=(AUTO/"p1b_vof_fullstress_stderr.log").read_text().splitlines()
auto_blocks=[]; current=None
for line in auto_stdout:
    if line.startswith("t mass_drift "):
        if current is not None: auto_blocks.append(current)
        current=[]
    elif current is not None:
        fields=line.split()
        try: vals=[float(x) for x in fields[:7]]
        except (ValueError, IndexError): continue
        if len(vals)==7: current.append(vals)
if current is not None: auto_blocks.append(current)
assert len(auto_blocks)==3, f"expected 3 automatic-curvature blocks, got {len(auto_blocks)}"
auto_mesh=[]
for cells, rows in zip((8,16,32),auto_blocks):
    m=next((re.match(rf"^{cells} ([0-9.eE+-]+) ([0-9.eE+-]+)$",s.strip()) for s in auto_stderr if re.match(rf"^{cells} ([0-9.eE+-]+) ([0-9.eE+-]+)$",s.strip())),None)
    assert m, f"missing automatic-curvature final speed for R/dx={cells}"
    measured_over_U=float(m.group(1)); theory_over_U=float(m.group(2)); ratio=measured_over_U/theory_over_U
    mass_pct=[100*r[1] for r in rows]
    auto_mesh.append({"cells_per_radius":cells,"N_domain":cells*16,"final_sample_t_over_t0":f"{rows[-1][0]:.9g}","measured_velocity_over_U":f"{measured_over_U:.9g}","theory_velocity_over_U":f"{theory_over_U:.9g}","measured_over_theory":f"{ratio:.9g}","relative_velocity_error_pct":f"{100*abs(ratio-1):.6f}","final_logged_mass_drift_pct":f"{mass_pct[-1]:.9g}","max_abs_logged_mass_drift_pct":f"{max(abs(x) for x in mass_pct):.9g}","log_samples":len(rows),"pass_speed_trend":""})
for i,row in enumerate(auto_mesh):
    row["pass_speed_trend"] = "coarsest baseline" if i==0 else ("error decreases with refinement" if float(row["relative_velocity_error_pct"]) < float(auto_mesh[i-1]["relative_velocity_error_pct"]) else "error did not decrease")
write_csv(AUTO/"p1b_autocurvature_velocity_grid_convergence.csv",auto_mesh)

# Combined deliverables make method-to-method differences explicit.
combined=[]; mass_rows=[]
for method, rows in (("official_CLSVOF_integral",mesh),("project_automatic_curvature_standard_VOF",auto_mesh)):
    for r in rows:
        row={"method":method,**r}; combined.append(row)
        mass_rows.append({"method":method,"cells_per_radius":r["cells_per_radius"],"final_logged_mass_drift_pct":r["final_logged_mass_drift_pct"],"max_abs_logged_mass_drift_pct":r["max_abs_logged_mass_drift_pct"],"log_samples":r["log_samples"]})
write_csv(HERE/"p1b_velocity_convergence.csv",combined)
write_csv(HERE/"p1b_mass_conservation.csv",mass_rows)

# User's supplied hand-extracted JSON is an explicitly non-original convenience subset.
summary={"doi":"10.1007/s10765-016-2104-7","raw_nist_doi_match":True,"raw_original_json_points_ag":len(raw_obs),"user_extracted_not_original_matches_raw_ag_subset":True,"temperature_range_K":[min(r[0] for r in raw_obs),max(r[0] for r in raw_obs)],"nist_property":"Viscosity, Pa*s","expanded_uncertainty_confidence_level_pct":95,"expanded_uncertainty_evaluation_method":"Propagation of evaluated standard uncertainties","fit":"ln(mu/(mPa*s)) = -0.4631 + 2089/T; mu from exp() in mPa*s","max_abs_standardized_fit_residual":max(abs(float(r["residual_over_expanded_uncertainty"])) for r in fit_rows),"max_abs_relative_fit_deviation_pct":max(abs(float(r["relative_fit_deviation_pct"])) for r in fit_rows),"warning":"The correlation is checked only over the NIST observation range; extrapolation at 2000 K and above is not experimentally validated by these 9 points."}
(NIST/"nist_ag_uncertainty_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps({"nist_points":len(raw_obs),"max_abs_fit_residual_in_U":summary["max_abs_standardized_fit_residual"],"max_abs_fit_deviation_pct":summary["max_abs_relative_fit_deviation_pct"],"p1b_official_clsvof":mesh,"p1b_project_autocurvature_vof":auto_mesh},indent=2))
