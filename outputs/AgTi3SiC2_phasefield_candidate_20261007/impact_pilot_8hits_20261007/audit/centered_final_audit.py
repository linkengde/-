from pathlib import Path
import json
import sys
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'generation'))
import map_atomic_candidate as ma

LX, LY, XY = 37.28597486102, 32.290601434511, -18.64298743051
CX, CY = 9.9573898076, 24.2578814220
SHIFTS = [(i * LX + j * XY, j * LY, 0.0) for i in (-1, 0, 1) for j in (-1, 0, 1)]
MASSES = np.array([47.867, 28.0855, 12.0107, 107.8682])
KB = 8.617333262e-5
AMU_A2PS2_TO_EV = 1.0364269e-4

def read_data(path):
    lines = Path(path).read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().startswith('Atoms')) + 1
    rows = []
    for line in lines[start:]:
        parts = line.split('#', 1)[0].split()
        if len(parts) < 5:
            if rows:
                break
            continue
        try:
            rows.append((int(parts[0]), int(parts[1]), *map(float, parts[2:5])))
        except ValueError:
            if rows:
                break
    a = np.asarray(rows, dtype=float)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0].astype(int), a[:, 1].astype(int), a[:, 2:5]

def read_last_dump(path):
    lines = Path(path).read_text().splitlines()
    frames = []
    k = 0
    while k < len(lines):
        if not lines[k].startswith('ITEM: TIMESTEP'):
            k += 1
            continue
        step = int(lines[k + 1])
        n = int(lines[k + 3])
        j = k + 9
        rows = np.asarray([[float(v) for v in lines[q].split()] for q in range(j, j + n)])
        frames.append((step, rows))
        k = j + n
    return frames[-1]

def radial(x):
    v = (x[:, 1] / LY) % 1.0
    u = ((x[:, 0] - XY * (x[:, 1] / LY)) / LX) % 1.0
    cv = CY / LY
    cu = ((CX - XY * cv) / LX) % 1.0
    du = (u - cu + 0.5) % 1.0 - 0.5
    dv = (v - cv + 0.5) % 1.0 - 0.5
    return np.sqrt((du * LX + dv * XY) ** 2 + (dv * LY) ** 2)

def pair_stats(ids_a, a, ids_b, b, same=False, radius=3.5):
    images = np.concatenate([b + np.asarray(s) for s in SHIFTS], axis=0)
    image_ids = np.tile(ids_b, len(SHIFTS))
    tree = cKDTree(images)
    dmin, idx = tree.query(a, k=1)
    min_i = int(np.argmin(dmin))
    result = {
        'min_A': float(dmin[min_i]),
        'pair_ids': [int(ids_a[min_i]), int(image_ids[idx[min_i]])],
        'below_1p5_A': 0,
        'below_2p0_A': 0,
        'below_2p5_A': 0,
    }
    pairs = {}
    for ia, point in enumerate(a):
        near = tree.query_ball_point(point, radius)
        for image_idx in near:
            ibid = int(image_ids[image_idx])
            aid = int(ids_a[ia])
            if same and aid == ibid:
                continue
            key = (min(aid, ibid), max(aid, ibid)) if same else (aid, ibid)
            dist = float(np.linalg.norm(point - images[image_idx]))
            if key not in pairs or dist < pairs[key]:
                pairs[key] = dist
    distances = np.fromiter(pairs.values(), dtype=float)
    if same and len(distances):
        best = min(pairs, key=pairs.get)
        result['min_A'] = float(pairs[best])
        result['pair_ids'] = [int(best[0]), int(best[1])]
    result['pairs_within_3p5_A'] = len(distances)
    for label, cutoff in [('below_1p5_A', 1.5), ('below_2p0_A', 2.0), ('below_2p5_A', 2.5)]:
        result[label] = int(np.sum(distances < cutoff))
    return result

def local_temperature(arr, rmax, only_ag=False):
    ids = arr[:, 0].astype(int)
    types = arr[:, 1].astype(int)
    x = arr[:, 2:5]
    v = arr[:, 5:8]
    mask = (ids <= 4305) & (x[:, 2] >= 12.0) & (radial(x) <= rmax)
    if only_ag:
        mask &= types == 4
    n = int(mask.sum())
    if n < 2:
        return {'n': n, 'T_K': None}
    mass = MASSES[types[mask] - 1]
    vv = v[mask] - np.average(v[mask], axis=0, weights=mass)
    ke = 0.5 * np.sum(mass[:, None] * vv ** 2) * AMU_A2PS2_TO_EV
    return {'n': n, 'T_K': float(2 * ke / ((3 * n - 3) * KB))}

def surface_profile(ids, x):
    r = radial(x)
    core = ids <= 4305
    surf = core & (x[:, 2] > 40.0)
    out = {}
    for lo, hi in ((0, 2), (2, 4), (4, 6), (6, 8), (8, 10)):
        z = x[surf & (r >= lo) & (r < hi), 2]
        out[f'{lo}-{hi}_A'] = {
            'n': int(len(z)),
            'max_z_A': float(z.max()) if len(z) else None,
            'p90_z_A': float(np.percentile(z, 90)) if len(z) else None,
            'median_z_A': float(np.median(z)) if len(z) else None,
        }
    out['center_n_z_gt40_A'] = int((core & (r <= 2) & (x[:, 2] > 40)).sum())
    out['center_n_z_gt45_A'] = int((core & (r <= 2) & (x[:, 2] > 45)).sum())
    return out

def centered_displacements(start, end):
    ids0, t0, x0 = start
    ids1, t1, x1 = end
    lookup = {int(i): j for j, i in enumerate(ids0)}
    keep = (ids1 <= 4305)
    match = np.array([lookup[int(i)] for i in ids1[keep]])
    delta = x1[keep] - x0[match]
    # Minimum image in the triclinic lateral cell; z remains non-periodic.
    dv = delta[:, 1] / LY
    du = (delta[:, 0] - XY * dv) / LX
    du -= np.rint(du)
    dv -= np.rint(dv)
    delta[:, 0] = du * LX + dv * XY
    delta[:, 1] = dv * LY
    distance = np.linalg.norm(delta, axis=1)
    rr = radial(x1[keep])
    types = t1[keep]
    out = {}
    for name, phase in [('TSC', types != 4), ('Ag_substrate', types == 4)]:
        local = phase & (rr <= 8.0)
        out[name] = {
            'n': int(phase.sum()),
            'all_rms_A': float(np.sqrt(np.mean(distance[phase] ** 2))),
            'all_p95_A': float(np.percentile(distance[phase], 95)),
            'all_gt1A': int((distance[phase] > 1.0).sum()),
            'r_le8_rms_A': float(np.sqrt(np.mean(distance[local] ** 2))),
            'r_le8_p95_A': float(np.percentile(distance[local], 95)),
            'r_le8_gt1A': int((distance[local] > 1.0).sum()),
            'r_le8_gt2A': int((distance[local] > 2.0).sum()),
        }
    return out

hit_path = ROOT / 'structures/impact_final_cool_0p10ps.data'
ctl_path = ROOT / 'structures/nohit_final_cool_0p10ps.data'
hit_ids, hit_t, hit_x = read_data(hit_path)
ctl_ids, ctl_t, ctl_x = read_data(ctl_path)
hit_step, hit_frame = read_last_dump(ROOT / 'dumps/impact_final_cool_0p10ps.dump')
ctl_step, ctl_frame = read_last_dump(ROOT / 'dumps/nohit_final_cool_0p10ps.dump')

hit_core = hit_ids <= 4305
ctl_core = ctl_ids <= 4305
hit_ag = hit_t == 4
ctl_ag = ctl_t == 4
hit_tsc = hit_t != 4
ctl_tsc = ctl_t != 4

# Parse the last LAMMPS thermo row, whose c_Tfree is the spatial free-zone estimate.
def final_tfree(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if len(p) >= 6 and p[0].isdigit():
            try:
                rows.append((int(p[0]), float(p[3]), int(p[2])))
            except ValueError:
                pass
    step, temp, atoms = rows[-1]
    return {'step': step, 'Tfree_K': temp, 'atoms': atoms}

out = {
    'model': {
        'substrate_atoms': 4305,
        'substrate_type_counts': {'Ti': 1464, 'Si': 519, 'C': 925, 'Ag': 1397},
        'cell_A': {'Lx': LX, 'Ly': LY, 'xy_tilt': XY, 'zlo': -10.0, 'zhi': 100.0},
        'boundaries': 'p p f',
        'actual_impact_center_A': [CX, CY],
        'thermal_region_center_A': [CX, CY],
    },
    'run': {'Ag_impacts': 4, 'energy_each_eV': 62.0, 'total_incident_eV': 248.0, 'total_time_ps': 0.94, 'last_cooling_window_ps': 0.10, 'final_dump_local_step': hit_step},
    'counts': {
        'hit_total': int(len(hit_ids)), 'hit_substrate': int(hit_core.sum()), 'hit_type_counts': {str(int(t)): int(np.sum(hit_t == t)) for t in np.unique(hit_t)},
        'control_total': int(len(ctl_ids)), 'control_type_counts': {str(int(t)): int(np.sum(ctl_t == t)) for t in np.unique(ctl_t)},
        'hit_ids_contiguous': bool(np.array_equal(np.sort(hit_ids), np.arange(1, len(hit_ids) + 1))),
        'control_ids_contiguous': bool(np.array_equal(np.sort(ctl_ids), np.arange(1, len(ctl_ids) + 1))),
    },
    'surface': {'hit': surface_profile(hit_ids, hit_x), 'control': surface_profile(ctl_ids, ctl_x)},
    'connectivity': {
        'cutoffs_A': {'TSC': 3.35, 'Ag': 3.2},
        'hit_TSC': ma.component_summary(hit_x[hit_tsc], 3.35, ma.cell_translations(LX, LY, XY)),
        'control_TSC': ma.component_summary(ctl_x[ctl_tsc], 3.35, ma.cell_translations(LX, LY, XY)),
        'hit_Ag_substrate': ma.component_summary(hit_x[hit_ag & hit_core], 3.2, ma.cell_translations(LX, LY, XY)),
        'control_Ag': ma.component_summary(ctl_x[ctl_ag], 3.2, ma.cell_translations(LX, LY, XY)),
    },
    'final_Tfree': {
        'hit': final_tfree(ROOT / 'logs/impact_final_cool_0p10ps.lammps'),
        'control': final_tfree(ROOT / 'logs/nohit_final_cool_0p10ps.lammps'),
    },
    'instantaneous_temperature_from_final_dump': {
        'hit_all_r4': local_temperature(hit_frame, 4.0), 'hit_all_r8': local_temperature(hit_frame, 8.0),
        'hit_Ag_r4': local_temperature(hit_frame, 4.0, True), 'hit_Ag_r8': local_temperature(hit_frame, 8.0, True),
        'control_all_r4': local_temperature(ctl_frame, 4.0), 'control_all_r8': local_temperature(ctl_frame, 8.0),
        'control_Ag_r4': local_temperature(ctl_frame, 4.0, True), 'control_Ag_r8': local_temperature(ctl_frame, 8.0, True),
    },
    'cross_phase_minima_hit': {
        f'{typ}-Ag': pair_stats(hit_ids[(hit_t == typ) & hit_core], hit_x[(hit_t == typ) & hit_core], hit_ids[hit_ag], hit_x[hit_ag])
        for typ in (1, 2, 3)
    },
    'Ag_Ag_min_hit': pair_stats(hit_ids[hit_ag], hit_x[hit_ag], hit_ids[hit_ag], hit_x[hit_ag], same=True),
    'min_Ti_C_hit': pair_stats(hit_ids[hit_t == 1], hit_x[hit_t == 1], hit_ids[hit_t == 3], hit_x[hit_t == 3]),
    'min_C_C_hit': pair_stats(hit_ids[hit_t == 3], hit_x[hit_t == 3], hit_ids[hit_t == 3], hit_x[hit_t == 3], same=True),
    'mobility_last_0p10ps': centered_displacements(read_data(ROOT / 'structures/impact_4th_Ag_62eV.data'), (hit_ids, hit_t, hit_x)),
    'dump_force_velocity_finite': bool(np.isfinite(hit_frame).all() and np.isfinite(ctl_frame).all()),
}

out_path = ROOT / 'audit/centered_final_audit.json'
out_path.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
