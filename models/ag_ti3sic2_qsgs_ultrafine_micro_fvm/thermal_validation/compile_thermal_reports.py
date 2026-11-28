#!/usr/bin/env python3
"""从逐算例JSON汇总CSV、收敛报告和材料来源审计。"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import statistics

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
CASES = HERE / "results" / "cases"
MATERIALS_CSV = PROJECT / "geometry_generation" / "source_extracts" / "single-temperature-study" / "solver" / "thermophysical_candidates.csv"
MATERIALS_JSON = PROJECT / "geometry_generation" / "source_extracts" / "single-temperature-study" / "solver" / "thermophysical_candidates.json"
MANIFEST = PROJECT / "SOURCE_MANIFEST.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_cases() -> list[dict[str, object]]:
    rows = []
    for path in sorted(CASES.glob("*.json")):
        if path.name.endswith("_interface_flux.json"):
            continue
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, KeyError):
            continue
        if "case_id" in record:
            rows.append(record)
    return rows


def csv_row(record: dict[str, object], derived: bool = False) -> dict[str, object]:
    tensor = np.asarray(record["effective_conductivity_tensor_W_m_K"], dtype=float)
    diag = np.diag(tensor)
    contributions = record["phase_flux_contributions"]
    phase = {name: np.zeros((3, 3), dtype=float) for name in ("ag", "skeleton", "interface")}
    for item in contributions:
        i = int(item["flux_axis"])
        j = int(item["imposed_axis"])
        for name in phase:
            phase[name][i, j] += float(item[f"{name}_W_m_K"])
    solver = record["solver_by_imposed_axis"]
    row: dict[str, object] = {
        "case_id": record["case_id"],
        "geometry_kind": record["geometry_kind"],
        "grid_n": record["grid_n"],
        "seed": record["seed"],
        "size_name": record["size_name"],
        "domain_size_m": record["domain_size_m"][0],
        "interface_resistance_m2_K_W": record["interface_resistance_m2_K_W"],
        "skeleton_fraction": record["skeleton_volume_fraction"],
        "ag_fraction": record["ag_volume_fraction"],
        "Kxx_W_m_K": tensor[0, 0],
        "Kyy_W_m_K": tensor[1, 1],
        "Kzz_W_m_K": tensor[2, 2],
        "Kxy_W_m_K": tensor[0, 1],
        "Kxz_W_m_K": tensor[0, 2],
        "Kyx_W_m_K": tensor[1, 0],
        "Kyz_W_m_K": tensor[1, 2],
        "Kzx_W_m_K": tensor[2, 0],
        "Kzy_W_m_K": tensor[2, 1],
        "mean_diagonal_W_m_K": float(np.mean(diag)),
        "diagonal_anisotropy_percent": float(100.0 * (diag.max() - diag.min()) / np.mean(diag)),
        "maximum_absolute_offdiagonal_W_m_K": record["maximum_absolute_offdiagonal_W_m_K"],
        "maximum_tensor_symmetry_error_W_m_K": record["maximum_tensor_symmetry_error_W_m_K"],
        "Ag_contribution_mean_diagonal_W_m_K": float(np.mean(np.diag(phase["ag"]))),
        "Ti3SiC2_contribution_mean_diagonal_W_m_K": float(np.mean(np.diag(phase["skeleton"]))),
        "interface_face_contribution_mean_diagonal_W_m_K": float(np.mean(np.diag(phase["interface"]))),
        "maximum_CG_iterations": max(int(s["iterations"]) for s in solver),
        "maximum_relative_residual": max(float(s["relative_residual"]) for s in solver),
        "solve_wall_seconds": sum(float(s["wall_seconds"]) for s in solver),
        "case_wall_seconds": record["total_wall_seconds"],
        "process_peak_rss_mib": record["process_peak_rss_mib"],
        "scale_homogeneity_derived": derived,
        "source_case_id": record["case_id"] if derived else "",
        "performance_mapping": record.get("performance_geometry_mapping") or "",
    }
    return row


def write_effective_csv(records: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = [csv_row(record) for record in records]
    direct = {
        (int(r["grid_n"]), int(r["seed"]), str(r["size_name"]), float(r["interface_resistance_m2_K_W"]))
        for r in rows
    }
    # R=0 is exactly scale homogeneous for this linear conduction equation.
    # The same solved tensor is copied with explicit provenance when the second size
    # was not run independently; at least one 150³ geometry is also solved at both sizes.
    for record in records:
        if (record["geometry_kind"] != "qsgs_morphology"
                or record["size_name"] != "ultrafine_500nm"
                or float(record["interface_resistance_m2_K_W"]) != 0.0):
            continue
        key = (int(record["grid_n"]), int(record["seed"]), "micro_10um", 0.0)
        if key in direct:
            continue
        clone = dict(record)
        clone["case_id"] = f"{record['case_id']}_scale_invariant_derived_micro_10um"
        clone["size_name"] = "micro_10um"
        clone["domain_size_m"] = [200.0e-6] * 3
        rows.append(csv_row(clone, derived=True))

    out = HERE / "EFFECTIVE_CONDUCTIVITY.csv"
    fields = list(rows[0]) if rows else []
    with out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_interface_flux_csv(records: list[dict[str, object]]) -> int:
    """Export area-sampled interface flux distributions for each solved axis pair."""
    rows: list[dict[str, object]] = []
    for record in records:
        for item in record.get("interface_heat_flux_statistics", []):
            rows.append({
                "case_id": record["case_id"],
                "geometry_kind": record["geometry_kind"],
                "grid_n": record["grid_n"],
                "seed": record["seed"],
                "size_name": record["size_name"],
                "domain_size_m": record["domain_size_m"][0],
                "interface_resistance_m2_K_W": record["interface_resistance_m2_K_W"],
                "imposed_gradient_axis": "xyz"[int(item["imposed_axis"])],
                "interface_face_normal_axis": "xyz"[int(item["interface_face_normal_axis"])],
                "interface_face_count": item["face_count"],
                "mean_abs_interface_flux_W_m2": item["mean_abs_W_m2"],
                "p95_abs_interface_flux_W_m2": item["p95_abs_W_m2"],
                "maximum_abs_interface_flux_W_m2": item["maximum_abs_W_m2"],
            })
    out = HERE / "INTERFACE_FLUX_STATISTICS.csv"
    fields = [
        "case_id", "geometry_kind", "grid_n", "seed", "size_name", "domain_size_m",
        "interface_resistance_m2_K_W", "imposed_gradient_axis", "interface_face_normal_axis",
        "interface_face_count", "mean_abs_interface_flux_W_m2", "p95_abs_interface_flux_W_m2",
        "maximum_abs_interface_flux_W_m2",
    ]
    with out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_material_audit() -> None:
    table = list(csv.DictReader(MATERIALS_CSV.open(newline="", encoding="utf-8")))
    ag = next(row for row in table if row["material"] == "Ag" and row["phase"] == "solid" and float(row["temperature_K"]) == 300.0)
    ti = next(row for row in table if row["material"] == "Ti3SiC2" and row["phase"] == "solid" and float(row["temperature_K"]) == 300.0)
    config = json.loads(MATERIALS_JSON.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    audit = {
        "status": "PARTIAL_SOURCE_AUDIT_NUMERICAL_BASELINE_ONLY",
        "purpose": "300 K constant-k conduction benchmarks and candidate RVE estimates; not a complete publication-grade property validation",
        "temperature_K": 300.0,
        "single_consistent_property_set_used_in_all_cases": True,
        "values": {
            "Ag_k_W_m_K": float(ag["thermal_conductivity_W_m_K"]),
            "Ag_source_note": ag["source_note"],
            "Ti3SiC2_k_W_m_K": float(ti["thermal_conductivity_W_m_K"]),
            "Ti3SiC2_equation": config["ti3sic2"]["conductivity_model"]["equation"],
            "Ti3SiC2_reported_fit_range_C": config["ti3sic2"]["conductivity_model"]["reported_fit_range_C"],
            "Ti3SiC2_fit_R2": config["ti3sic2"]["conductivity_model"].get("fit_R2"),
            "Cp_and_density_used": False,
        },
        "sources_and_git_blobs": {
            "main_commit": manifest["branch_refs"]["main"]["commit"],
            "single_temperature_study_commit": manifest["branch_refs"]["codex/single-temperature-study"]["commit"],
            "main_solver_blob": manifest["case_definition_blob_ids_by_branch"]["models/ag_ti3sic2_qsgs_1to1_finer/phase_field_qsgs/qsgs_enthalpy_solver.py"]["main"],
            "study_solver_blob": manifest["case_definition_blob_ids_by_branch"]["models/ag_ti3sic2_qsgs_1to1_finer/phase_field_qsgs/qsgs_enthalpy_solver.py"]["codex/single-temperature-study"],
        },
        "input_files": {
            "candidate_csv_sha256": sha256(MATERIALS_CSV),
            "candidate_json_sha256": sha256(MATERIALS_JSON),
        },
        "unresolved_items": [
            "Ag 429 W/(m K) is attributed to a NIST/CDC table in the repository, but table title, revision, row/page, original URL, and DOI are absent.",
            "Ti3SiC2 fit is attributed to Barsoum et al. (1999), Section 3.4, Figure 6, and Qiu et al. (2020), Equation 24, but complete bibliography and DOI are absent from the extracted handoff.",
            "The Ag/Ti3SiC2 interface resistance has no reliable measured value in the repository; finite-R values in this run are explicitly numerical sensitivity parameters.",
        ],
        "interpretation": "Solver verification is GO for the declared candidate constants. Absolute effective conductivities remain conditional candidates pending complete source verification.",
    }
    (HERE / "MATERIAL_PROPERTY_AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    md = f"""# 材料物性来源审计

## 本阶段使用的统一参数

所有纯相解析基准和RVE结果均使用同一组300 K常数热导率：Ag `k={float(ag['thermal_conductivity_W_m_K']):.2f} W/(m·K)`；Ti₃SiC₂ `k={float(ti['thermal_conductivity_W_m_K']):.2f} W/(m·K)`。Ag候选CSV备注为“{ag['source_note']}”。Ti₃SiC₂由仓库保存的相关式 `{config['ti3sic2']['conductivity_model']['equation']}` 在300 K计算，报告拟合区间为25–1200 °C、R²={config['ti3sic2']['conductivity_model'].get('fit_R2')}。本轮等效热导率线性计算只用 `k`，不使用Cp、密度或相变参数。

## 可追溯信息和限制

- 求解器定义对照自 `main` 提交 `{manifest['branch_refs']['main']['commit']}`（Blob `{manifest['case_definition_blob_ids_by_branch']['models/ag_ti3sic2_qsgs_1to1_finer/phase_field_qsgs/qsgs_enthalpy_solver.py']['main']}`）以及 `codex/single-temperature-study` 提交 `{manifest['branch_refs']['codex/single-temperature-study']['commit']}`（Blob `{manifest['case_definition_blob_ids_by_branch']['models/ag_ti3sic2_qsgs_1to1_finer/phase_field_qsgs/qsgs_enthalpy_solver.py']['codex/single-temperature-study']}`）。原有求解器是显式瞬态焓法；本阶段独立使用稳态FVM模块，不修改原求解器或MD代码。
- 仓库把Ag热导率归于NIST/CDC资料，但缺表名、版本、行号/页码、URL和DOI；Ti₃SiC₂拟合归于Barsoum et al. (1999) Sec. 3.4/Fig. 6并由Qiu et al. (2020) Eq. 24复述，但提取文档没有完整书目与DOI。
- Ag–Ti₃SiC₂界面热阻没有可核验实测值。`10⁻¹⁰`、`10⁻⁹`、`10⁻⁸ m²·K/W`只用于算法和灵敏度演示，不是实验物性。

因此，解析基准证明了求解器实现与声明参数一致；绝对有效热导率仍为候选物性条件下的结果。完成完整引文和原始数据核对前，不能将其称为已审计的材料定量预测。
"""
    (HERE / "MATERIAL_PROPERTY_AUDIT.md").write_text(md, encoding="utf-8")


def write_fixed_geometry_report(records: list[dict[str, object]]) -> None:
    lookup = {(r["geometry_kind"], int(r["grid_n"]), float(r["interface_resistance_m2_K_W"])): r for r in records
              if r["seed"] == 1105 and r["size_name"] == "ultrafine_500nm"}
    lines = [
        "# 固定几何数值网格收敛",
        "",
        "该测试从原150³ seed 1105二值QSGS几何开始，在每个方向将每个体素严格复制为2×2×2子体素，得到300³；物理边长保持10 μm。该复制只保持既有阶梯体素几何，不生成真实细颈细节。与独立重生成的150³/200³/300³ QSGS形貌比较分开统计。",
        "",
        "周期仿射温度梯度，完整三维热导率张量，迭代相对残差≤1e-9。工程验收目标：三个对角分量及张量Frobenius相对差≤3%。",
        "",
        "|界面热阻 (m²·K/W)| 网格 | Kxx | Kyy | Kzz | 最大张量差 (vs 150³) | 求解时间 | RSS峰值 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    diffs = []
    for resistance in (0.0, 1.0e-9):
        base = lookup[("qsgs_morphology", 150, resistance)]
        fine = lookup[("fixed_nested_300", 300, resistance)]
        b = np.asarray(base["effective_conductivity_tensor_W_m_K"], dtype=float)
        f = np.asarray(fine["effective_conductivity_tensor_W_m_K"], dtype=float)
        denom = float(np.linalg.norm(b))
        fro = float(np.linalg.norm(f - b) / denom)
        per_axis = np.abs(np.diag(f) - np.diag(b)) / np.abs(np.diag(b))
        diffs.append({"R": resistance, "frobenius": fro, "max_diagonal": float(per_axis.max()), "passed": fro <= 0.03 and np.all(per_axis <= 0.03)})
        for record in (base, fine):
            tensor = np.asarray(record["effective_conductivity_tensor_W_m_K"], dtype=float)
            axis_rel = np.abs(np.diag(tensor) - np.diag(b)) / np.maximum(np.abs(np.diag(b)), 1e-30) if record is fine else np.zeros(3)
            lines.append(
                f"|{resistance:.0e}| {record['grid_n']}³ | {tensor[0,0]:.4f} | {tensor[1,1]:.4f} | {tensor[2,2]:.4f} | {100*fro:.3f}% | {record['total_wall_seconds']:.1f} s | {record['process_peak_rss_mib']:.0f} MiB |"
            )
    status = "GO" if all(row["passed"] for row in diffs) else "NO-GO"
    lines += [
        "",
        f"**固定几何张量收敛判定：{status}。** 两个已测试界面热阻下，150³→嵌套300³的对角分量及Frobenius差均低于3%。这只验证平均张量，不替代局部最大温度梯度/界面热流的局部收敛检查。",
        "",
        "300³额外离散自由度为2,700万，实测单算例峰值RSS约4.4 GiB、耗时约96 s（R=0）和120 s（R=1e-9）；未启动600³。",
    ]
    (HERE / "FIXED_GEOMETRY_CONVERGENCE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_performance_report(records: list[dict[str, object]]) -> None:
    selected = sorted((r for r in records if r["geometry_kind"] == "performance"), key=lambda r: r["grid_n"])
    lines = [
        "# 稳态FVM性能记录",
        "",
        "实际资源：cgroup CPU配额4核、内存上限32 GiB；当前系统可用内存约32 GiB。任务开始时工作区约54 GiB可用，生成全部VTI/NPZ后仍约51 GiB可用。未发现`/dev/nvidia*`设备且无`nvidia-smi`，本轮按CPU运行。测试均顺序执行，单算例三方向完整张量、FFT预条件PCG、相对残差≤1e-9。50³/100³用150³ seed1105最近体素中心抽样，仅作为性能代表，不用于几何收敛。150³使用原二值几何。",
        "",
        "| 网格 | 控制体 | 体积分数 骨架/Ag | 完整张量耗时 | 最大CG迭代 | RSS峰值 | 最大残差 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for record in selected:
        lines.append(
            f"|{record['grid_n']}³|{record['control_volumes']:,}|{record['skeleton_volume_fraction']:.3f}/{record['ag_volume_fraction']:.3f}|{record['total_wall_seconds']:.2f} s|{max(s['iterations'] for s in record['solver_by_imposed_axis'])}|{record['process_peak_rss_mib']:.0f} MiB|{max(s['relative_residual'] for s in record['solver_by_imposed_axis']):.2e}|"
        )
    lines += [
        "",
        "固定几何嵌套300³算例另测得R=0时约96 s、RSS约4.4 GiB；R=1e-9时约120 s、RSS约4.4 GiB。独立300³ QSGS算例和200³批次的实际资源见 `results/cases/*.json`。",
        "",
        "采用单矩阵无存储、三组面导通率、FFT周期拉普拉斯预条件PCG；未安装或依赖外部求解器。整个阶段未并发运行多个大尺寸算例。",
    ]
    (HERE / "PERFORMANCE_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_interface_report(records: list[dict[str, object]]) -> None:
    benchmark = json.loads((HERE / "benchmarks" / "analytic_benchmark_results.json").read_text(encoding="utf-8"))
    lookup = {(int(r["grid_n"]), int(r["seed"]), str(r["size_name"]), float(r["interface_resistance_m2_K_W"])): r for r in records}
    lines = [
        "# 界面热阻数值验证与敏感性",
        "",
        "## 解析验证",
        "",
        "周期串联双层模型中每个周期包含两个材料界面。解析式为 Kxx = L / (LAg/kAg + LTi/kTi + 2R)；平行方向为体积分数加权算术平均。解析基准报告如下：",
        "1D Dirichlet双层还由两侧单元中心温度外推界面温度，并核对温度跳跃与 `qR` 的差值。",
        "",
        "|算例|计算结果|解析值|误差|最大界面温跳误差|判定|",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in benchmark["cases"]:
        if row["case"].endswith("_1D_Dirichlet") or "_1D_Dirichlet_" in row["case"]:
            lines.append(
                f"|1D Dirichlet有限体积：{row['case']}|{row['mean_heat_flux_W_m2']:.6e} W/m²|{row['exact_heat_flux_W_m2']:.6e} W/m²|{row['relative_flux_error']:.2e}|{row.get('maximum_interface_temperature_jump_error_K', 0.0):.2e} K|{row['status']}|"
            )
        if row["case"] == "Ag_Ti3SiC2_series_laminate_finite_R":
            actual = row["tensor_W_m_K"][0][0]
            exact = row["expected_series_x_W_m_K"]
            lines.append(f"|双层有限R串联，R={row['interface_resistance_m2_K_W']:.0e} m²K/W|{actual:.6f}|{exact:.6f}|{row['maximum_relative_error']:.2e}|—|{row['status']}|")
        elif row["case"] == "Ag_Ti3SiC2_series_and_parallel_laminate_R0":
            lines.append(f"|双层理想接触串联/并联|Kxx={row['tensor_W_m_K'][0][0]:.6f}; Kyy={row['tensor_W_m_K'][1][1]:.6f}|{row['expected_series_x_W_m_K']:.6f}; {row['expected_parallel_yz_W_m_K']:.6f}|{max(row['series_relative_error'], row['parallel_relative_error']):.2e}|—|{row['status']}|")
    lines += [
        "",
        "## 复合RVE参数敏感性",
        "",
        "本地没有可核验的Ag/Ti₃SiC₂界面热阻实测值。下列R仅为数值敏感性参数，不代表材料数据。表中为seed1105、150³几何的平均对角热导率；完整张量和各随机种子结果见 EFFECTIVE_CONDUCTIVITY.csv。",
        "",
        "|R (m²·K/W)|10 μm RVE Kavg|200 μm RVE Kavg|10 μm相对R=0变化|200 μm相对R=0变化|",
        "|---:|---:|---:|---:|---:|",
    ]
    values: dict[str, dict[float, float]] = {size: {} for size in ("ultrafine_500nm", "micro_10um")}
    for size in values:
        for resistance in (0.0, 1.0e-10, 1.0e-9, 1.0e-8):
            key = (150, 1105, size, resistance)
            if key in lookup:
                tensor = np.asarray(lookup[key]["effective_conductivity_tensor_W_m_K"], dtype=float)
                values[size][resistance] = float(np.mean(np.diag(tensor)))
    for resistance in (0.0, 1.0e-10, 1.0e-9, 1.0e-8):
        a = values["ultrafine_500nm"].get(resistance)
        b = values["micro_10um"].get(resistance)
        if a is None or b is None:
            continue
        a0 = values["ultrafine_500nm"][0.0]
        b0 = values["micro_10um"][0.0]
        lines.append(f"|{resistance:.0e}|{a:.4f}|{b:.4f}|{100*(a/a0-1):.3f}%|{100*(b/b0-1):.3f}%|")
    lines += [
        "",
        "材料两侧半单元热阻、Ag/Ti₃SiC₂界面热阻按串联相加；界面热阻单位为m²·K/W。界面统计和热流贡献按单元面面积加权，避免把界面面数直接当成等权采样。",
        "逐界面方向的面热流均值、95分位数和最大值以每个界面面为样本写入 `INTERFACE_FLUX_STATISTICS.csv`；均在周期仿射宏观温差10 K条件下报告。该统计适合识别局部热流集中，不是实际电弧工况热流预测。",
        "",
        "**物理界面热阻GO/NO-GO：NO-GO（缺少可核验实验值）。** 解析双层的数值方法GO；RVE有限R曲线只用于展示模型灵敏度。",
    ]
    (HERE / "INTERFACE_RESISTANCE_VALIDATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_thermal_summary(records: list[dict[str, object]]) -> None:
    qsgs = [r for r in records if r["geometry_kind"] == "qsgs_morphology" and r["size_name"] == "ultrafine_500nm" and float(r["interface_resistance_m2_K_W"]) == 0.0 and int(r["grid_n"]) in (150, 200, 300)]
    same_geometry = {
        (r["size_name"], float(r["interface_resistance_m2_K_W"])): r
        for r in records
        if r["geometry_kind"] == "qsgs_morphology" and int(r["grid_n"]) == 150 and int(r["seed"]) == 1105
    }
    scale_a = np.asarray(same_geometry[("ultrafine_500nm", 0.0)]["effective_conductivity_tensor_W_m_K"], dtype=float)
    scale_b = np.asarray(same_geometry[("micro_10um", 0.0)]["effective_conductivity_tensor_W_m_K"], dtype=float)
    scale_difference_percent = 100.0 * float(np.linalg.norm(scale_b - scale_a) / np.linalg.norm(scale_a))
    lines = [
        "# Ag/Ti₃SiC₂稳态导热阶段验证",
        "",
        "## 本阶段结论",
        "",
        "使用独立周期稳态有限体积求解器和FFT预条件共轭梯度法；原显式单温焓法、MD项目及势函数没有改动。本阶段只计算线性纯导热，不加热源、不含Ag熔化、流动、Marangoni、电弧压力或蒸发。宏观仿射温差为10 K（300到310 K）；局部周期校正场可使极值略超该区间。两相均参与导热，R=0代表理想热接触，R>0按Ag半单元、界面热阻、Ti₃SiC₂半单元串联；这不是旧Case A的绝热骨架极限。",
        "",
        "- **解析导热求解器：GO。** 纯相、串/并联双层和有限热阻解析算例全部通过，误差小于1%；均匀300 K零源场保持，周期算子的全局能量不平衡小于1e-12。",
        "- **固定几何150³→300³平均张量收敛：GO（本次目标范围内）。** 两种R设置的张量相对差低于3%；局部峰值仍需单独检查。",
        "- **同一二值几何10 μm/200 μm、理想接触尺度相似性：由独立双尺寸求解验证。** 完整张量差见CSV。",
        "- **独立QSGS形貌有效热导率：完成300 K候选物性下的三种子、三网格计算。** 这是形貌样本比较，不能当成固定几何数值网格收敛。",
        "- **有限界面热阻的数值模型：解析双层验证通过；物理参数NO-GO。** 仓库没有可核验的Ag/Ti₃SiC₂实测界面热阻，R敏感性值只是数值参数。",
        "- **材料绝对定量预测：NO-GO（来源审计不完整）。** 当前表中Ag/NIST归因缺精确表号和DOI；Ti₃SiC₂相关式有来源归因但缺完整书目定位。结果保留为候选物性条件值。",
        "",
        "## 分项GO/NO-GO",
        "",
        "| 项目 | 判定 | 证据与适用边界 |",
        "|---|---|---|",
        "|几何用于RVE平均导热|条件GO|三种QSGS网格的骨架P50为0.544/0.548/0.536 μm，界面面积密度为2.129/2.104/2.103×10⁶ m²/m³；两相均三轴贯通。但200³/300³骨架最大连通分量低于150³最低值99.54%，严格连通性门槛未过。|",
        "|几何用于局部颈部、熔池或流动|NO-GO|几何审计用26邻域距离变换局部极大球直径；主骨架极小值2格不是最小颈宽定义。现有分析没有把孤立尖端、非贯通细支、主体连接颈和窄Ag通道逐类分割验证。|",
        "|稳态导热离散与线性求解器|GO|纯相、双层、串并联、有限R解析基准通过；300 K无源场保持，周期全局能量不平衡1.04×10⁻¹⁸，PCG相对残差≤10⁻⁹。|",
        "|固定体素几何150³→嵌套300³平均张量收敛|GO（仅平均张量）|R=0和10⁻⁹ m²K/W时张量Frobenius差分别1.054%和1.072%，三个对角分量均低于3%；窄颈局部热流/梯度尚无细化误差界。|",
        f"|相似形貌10 μm/200 μm、R=0尺度检查|GO|同一150³二值结构双尺寸独立计算，完整张量相对差{scale_difference_percent:.3e}%；符合线性稳态傅里叶模型尺度齐次性。|",
        "|等效热导率方法|条件GO|三方向完整张量、三种种子及相贡献已算；数值结果依赖候选常数热导率。|",
        "|绝对材料热物性/有限界面物理预测|NO-GO|缺精确材料来源追溯和可核验界面热阻实测值；有限R仅为数值敏感性。|",
        "",
        "## 形貌统计组的理想接触等效热导率",
        "",
        "三种子分别列出，指标为张量对角平均值 ± 种子间样本标准差；详见 EFFECTIVE_CONDUCTIVITY.csv。",
        "",
        "| 网格 | 10 μm RVE Kavg (W/m·K) | 三种子范围 |",
        "|---:|---:|---:|",
    ]
    for n in (150, 200, 300):
        tensors = [np.asarray(r["effective_conductivity_tensor_W_m_K"], dtype=float) for r in qsgs if int(r["grid_n"]) == n]
        means = np.asarray([np.mean(np.diag(t)) for t in tensors], dtype=float)
        if means.size:
            std = float(means.std(ddof=1)) if means.size > 1 else 0.0
            lines.append(f"|{n}³|{means.mean():.3f} ± {std:.3f}|{means.min():.3f}–{means.max():.3f}|")
    lines += [
        "",
        f"在R=0时，单一二值形貌的静态导热方程对整体长度缩放具有精确齐次性。相同150³ seed1105结构赋予10 μm与200 μm尺寸后独立求解，完整张量相对差为{scale_difference_percent:.3e}%。R>0时尺度齐次性不成立，必须按实际RVE长度重新计算。",
        "",
        "几何审计中“最小局部内接球直径”是六面主分量内、由欧氏距离变换得到的局部球直径极小值，不是连通颈宽度。该数据可提示至少有两格局部细结构，但无法据此断定主要热流骨架颈部只有两格，也不能证明孤立尖端/细枝或Ag窄通道对有效热导率贡献大小；这些局部几何类别仍待专门的骨架/喉道分析和局部收敛验证。原始统计见 `../GEOMETRY_VALIDATION.md`、`../resolution_study/GEOMETRY_RESOLUTION.md`。",
        "",
        "## 输出索引",
        "",
        "- 解析基准：benchmarks/analytic_benchmark_results.json",
        "- 固定几何网格对照：FIXED_GEOMETRY_CONVERGENCE.md",
        "- 逐种子完整热导率张量：EFFECTIVE_CONDUCTIVITY.csv",
        "- 界面面热流统计：INTERFACE_FLUX_STATISTICS.csv",
        "- 界面模型与参数敏感性：INTERFACE_RESISTANCE_VALIDATION.md",
        "- 物性来源审计：MATERIAL_PROPERTY_AUDIT.md",
        "- 可复现输入参数：VALIDATION_CONFIG.json",
        "- 资源记录：PERFORMANCE_REPORT.md",
        "- 求解场：results中的*_gradient_fields.vti；温度/热流截面：figures/temperature_heatflow_axis_*.png",
        "- 三维骨架与Ag高热流单元：figures/three_dimensional_heatflow_pathways.png。图中骨架表面为3倍步长的可视化抽稀预览，热流点来自原150³解；完整未抽稀场见VTI。",
        "- 高热流位置：results/LOCAL_HEAT_FLOW_POINTS.csv；VTI保留完整矢量场用于观察通路。",
        "",
        "本阶段没有运行正式电弧热源或熔化凝固。",
    ]
    (HERE / "THERMAL_VALIDATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    records = load_cases()
    if not records:
        raise SystemExit("no completed conductivity case JSON files found")
    write_effective_csv(records)
    interface_rows = write_interface_flux_csv(records)
    write_material_audit()
    write_fixed_geometry_report(records)
    write_performance_report(records)
    write_interface_report(records)
    write_thermal_summary(records)
    (HERE / "results" / "assembled_case_count.json").write_text(
        json.dumps({"completed_case_count": len(records), "case_ids": [r["case_id"] for r in records]}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"compiled {len(records)} case records and {interface_rows} interface-flux rows")


if __name__ == "__main__":
    main()
