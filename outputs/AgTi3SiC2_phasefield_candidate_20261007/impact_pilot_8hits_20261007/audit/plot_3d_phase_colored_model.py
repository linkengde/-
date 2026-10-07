from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LX, LY, XY = 37.28597486102, 32.290601434511, -18.64298743051
CX, CY = 9.9573898076, 24.2578814220
CV, CU = CY / LY, ((CX - XY * (CY / LY)) / LX) % 1


def read_data(path):
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().startswith('Atoms')) + 1
    rows = []
    for line in lines[start:]:
        p = line.split('#', 1)[0].split()
        if len(p) < 5:
            if rows: break
            continue
        try: rows.append((int(p[0]), int(p[1]), *map(float, p[2:5])))
        except ValueError:
            if rows: break
    a = np.asarray(rows)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0].astype(int), a[:, 1].astype(int), a[:, 2:5]


def local_xyz(xyz):
    v = (xyz[:, 1] / LY) % 1
    u = ((xyz[:, 0] - XY * (xyz[:, 1] / LY)) / LX) % 1
    du = (u - CU + 0.5) % 1 - 0.5
    dv = (v - CV + 0.5) % 1 - 0.5
    return np.column_stack((du * LX + dv * XY, dv * LY, xyz[:, 2]))

path = ROOT / 'structures/impact_post_cool_0p05ps.data'
ids, types, xyz = read_data(path)
coords = local_xyz(xyz)
phase = np.where(types == 4, 1, 0)  # 0 = Ti3SiC2, 1 = Ag
colors = {0: '#7850aa', 1: '#e4ba58'}
labels = {0: 'Ti₃SiC₂ skeleton (Ti, Si, C)', 1: 'Ag filling phase'}

fig = plt.figure(figsize=(13.5, 6.4), constrained_layout=True)
for panel, mode in enumerate(('full', 'cutaway'), 1):
    ax = fig.add_subplot(1, 2, panel, projection='3d')
    if mode == 'full':
        keep = np.ones(len(ids), dtype=bool)
        title = 'Full atomistic model'
    else:
        # Remove one upper-right quarter only for this visualization, exposing
        # the interior network without changing or editing the source model.
        keep = ~((coords[:, 0] > 1.5) & (coords[:, 1] > 1.5) & (ids <= 4305))
        title = 'Quarter cutaway: interior network exposed'
    for ph in (0, 1):
        m = keep & (phase == ph) & (ids <= 4305)
        ax.scatter(coords[m, 0], coords[m, 1], coords[m, 2],
                   s=5.0 if ph == 0 else 7.0, c=colors[ph], alpha=0.83,
                   depthshade=False, linewidths=0, label=labels[ph])
    m = keep & (ids > 4305)
    if m.any():
        ax.scatter(coords[m, 0], coords[m, 1], coords[m, 2], s=28,
                   c='#e33b33', marker='D', edgecolors='black', linewidths=.3,
                   depthshade=False, label='Deposited Ag atoms')
    ax.scatter([0], [0], [52], marker='x', s=55, color='red', linewidths=1.8)
    ax.set_xlim(-28, 28); ax.set_ylim(-17, 17); ax.set_zlim(-10, 54)
    ax.set_box_aspect((56, 34, 64))
    ax.view_init(elev=20, azim=-55)
    ax.set_xlabel('Local x (Å)', labelpad=3)
    ax.set_ylabel('Local y (Å)', labelpad=3)
    ax.set_zlabel('z (Å)', labelpad=3)
    ax.set_title(title, pad=8)
    ax.grid(False)
    if panel == 1: ax.legend(loc='upper left', fontsize=8, framealpha=.9)
fig.suptitle('The 3D crater model with both phases visible\nPurple = Ti₃SiC₂ skeleton; gold = Ag filling phase', fontsize=14)
out = ROOT / '3d_phase_colored_crater_model_8hits_1p85ps.png'
fig.savefig(out, dpi=200, bbox_inches='tight')
print(f'{out} {out.stat().st_size} bytes')
