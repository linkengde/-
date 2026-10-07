from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

LX, LY, XY = 37.28597486102, 32.290601434511, -18.64298743051
CX, CY = 9.9573898076, 24.2578814220
CV = CY / LY
CU = ((CX - XY * CV) / LX) % 1
ROOT = Path(__file__).resolve().parents[1]
SHIFTS = np.array([(i * LX + j * XY, j * LY, 0.0) for i in (-1, 0, 1) for j in (-1, 0, 1)])

def read_data(path):
    lines = Path(path).read_text().splitlines()
    i = next(j for j, s in enumerate(lines) if s.strip().startswith('Atoms')) + 1
    rows = []
    for line in lines[i:]:
        p = line.split('#', 1)[0].split()
        if len(p) < 5:
            if rows: break
            continue
        try: rows.append((int(p[0]), int(p[1]), float(p[2]), float(p[3]), float(p[4])))
        except ValueError:
            if rows: break
    a = np.asarray(rows)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0].astype(int), a[:, 1].astype(int), a[:, 2:5]

def local_xy(x):
    v = (x[:, 1] / LY) % 1
    u = ((x[:, 0] - XY * (x[:, 1] / LY)) / LX) % 1
    du = (u - CU + .5) % 1 - .5
    dv = (v - CV + .5) % 1 - .5
    return du * LX + dv * XY, dv * LY

colors = {1: '#7b4ab5', 2: '#4f9ca9', 3: '#333333', 4: '#e0b34d'}
labels = {1: 'Ti', 2: 'Si', 3: 'C', 4: 'Ag'}
paths = {
    'No-hit control (1.85 ps)': ROOT / 'structures/nohit_post_cool_0p05ps.data',
    'Eight 62 eV Ag hits (1.85 ps)': ROOT / 'structures/impact_post_cool_0p05ps.data',
}
fig = plt.figure(figsize=(15, 5.2), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05])
records = {}
for col, (title, path) in enumerate(paths.items()):
    ids, types, xyz = read_data(path)
    dx, dy = local_xy(xyz)
    sel = (np.abs(dy) < 2.5) & (np.abs(dx) < 10) & (xyz[:, 2] > 30)
    ax = fig.add_subplot(gs[0, col])
    for typ in (1, 2, 3, 4):
        m = sel & (types == typ) & (ids <= 4305)
        if m.any(): ax.scatter(dx[m], xyz[m, 2], s=18, c=colors[typ], label=labels[typ], alpha=.84, linewidths=0)
    m = sel & (ids > 4305)
    if m.any(): ax.scatter(dx[m], xyz[m, 2], s=34, c='#d62728', marker='D', label='Ag impactors', edgecolors='black', linewidths=.3)
    ax.axvline(0, color='#777777', lw=.7, ls=':')
    ax.set_xlim(-10, 10); ax.set_ylim(30, 54)
    ax.set_xlabel('Lateral distance from impact center (Å)')
    if col == 0: ax.set_ylabel('Height z (Å)')
    ax.set_title(title); ax.grid(alpha=.15)
    if col == 0: ax.legend(fontsize=8, loc='lower left', ncol=2, frameon=True)
    core = (ids <= 4305) & (np.hypot(dx, dy) <= 2)
    records['control' if col == 0 else 'impact'] = {
        'atoms': int(len(ids)), 'center_r_le_2A_top_z_A': float(xyz[core, 2].max()),
        'center_r_le_2A_substrate_atoms_z_gt_35A': int(np.sum(core & (xyz[:, 2] > 35))),
    }

ids, types, xyz = read_data(paths['Eight 62 eV Ag hits (1.85 ps)'])
dx, dy = local_xy(xyz)
m = (xyz[:, 2] > 35) & (ids <= 4305)
ax = fig.add_subplot(gs[0, 2])
sc = ax.scatter(dx[m], dy[m], c=xyz[m, 2], s=21, cmap='viridis', vmin=38, vmax=51, linewidths=0)
for r in (2, 4, 8): ax.add_patch(plt.Circle((0, 0), r, fill=False, color='white' if r < 5 else 'black', lw=.9, ls='--', alpha=.85))
ax.scatter([0], [0], s=85, marker='x', color='red', linewidths=2, label='impact center')
ax.set_aspect('equal'); ax.set_xlim(-10, 10); ax.set_ylim(-10, 10)
ax.set_xlabel('In-plane distance x (Å)'); ax.set_ylabel('In-plane distance y (Å)')
ax.set_title('Eight-hit top view\n(substrate atoms z > 35 Å)'); ax.legend(fontsize=8, loc='upper right')
fig.colorbar(sc, ax=ax, label='Height z (Å)', shrink=.78)
fig.suptitle('Ag–Ti₃SiC₂: crater-like surface after additional 0.05 ps cooling', fontsize=14)
out = ROOT / 'audit/surface_morphology_8hits_1p85ps.png'
fig.savefig(out, dpi=180, bbox_inches='tight')

# Nearest Ag-Ag contact with full in-plane periodic images.
ag_mask = types == 4
ag_ids, ag_xyz = ids[ag_mask], xyz[ag_mask]
images = np.concatenate([ag_xyz + shift for shift in SHIFTS], axis=0)
image_ids = np.tile(ag_ids, len(SHIFTS))
tree = cKDTree(images)
min_pair = {'distance_A': float('inf')}
for idx, point in enumerate(ag_xyz[:, :2]):
    d, j = tree.query(ag_xyz[idx], k=min(12, len(image_ids)))
    d, j = np.atleast_1d(d), np.atleast_1d(j)
    candidates = [(float(dd), int(image_ids[jj])) for dd, jj in zip(d, j) if image_ids[jj] != ag_ids[idx]]
    if candidates:
        dd, other = min(candidates)
        if dd < min_pair['distance_A']:
            min_pair = {'distance_A': dd, 'pair_ids': [int(ag_ids[idx]), other]}
records['impact']['minimum_Ag_Ag_distance_A'] = min_pair
records['image'] = {'file': out.name, 'width_px': int(fig.get_size_inches()[0] * fig.dpi), 'height_px': int(fig.get_size_inches()[1] * fig.dpi), 'bytes': int(out.stat().st_size)}
(ROOT / 'audit/postimpact_cooling_0p05ps_audit.json').write_text(json.dumps(records, indent=2) + '\n')
print(json.dumps(records, indent=2))
