from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import distance_transform_edt, gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
LX, LY, XY = 37.28597486102, 32.290601434511, -18.64298743051
CX, CY = 9.9573898076, 24.2578814220
CV, CU = CY / LY, ((CX - XY * (CY / LY)) / LX) % 1
BIN = 1.5
N = int(24 / BIN)


def read_data(path):
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().startswith('Atoms')) + 1
    rows = []
    for line in lines[start:]:
        fields = line.split('#', 1)[0].split()
        if len(fields) < 5:
            if rows: break
            continue
        try: rows.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
        except ValueError:
            if rows: break
    a = np.asarray(rows)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0].astype(int), a[:, 1].astype(int), a[:, 2:5]


def local_xy(xyz):
    v = (xyz[:, 1] / LY) % 1
    u = ((xyz[:, 0] - XY * (xyz[:, 1] / LY)) / LX) % 1
    du = (u - CU + 0.5) % 1 - 0.5
    dv = (v - CV + 0.5) % 1 - 0.5
    return du * LX + dv * XY, dv * LY


def surface_envelope(ids, xyz):
    dx, dy = local_xy(xyz)
    radius = np.hypot(dx, dy)
    use = (ids <= 4305) & (radius <= 12.0) & (xyz[:, 2] > 28.0)
    ix = np.floor((dx[use] + 12.0) / BIN).astype(int)
    iy = np.floor((dy[use] + 12.0) / BIN).astype(int)
    source = np.flatnonzero(use)
    height = np.full((N, N), np.nan)
    atoms = {}
    for idx, i, j in zip(source, ix, iy):
        if not (0 <= i < N and 0 <= j < N):
            continue
        if np.isnan(height[j, i]) or xyz[idx, 2] > height[j, i]:
            height[j, i] = xyz[idx, 2]
            atoms[(j, i)] = int(idx)
    valid = np.isfinite(height)
    # Fill empty cells from the nearest sampled surface atom, then smooth only
    # the visualization envelope. The colored dots remain actual atom sites.
    nearest = distance_transform_edt(~valid, return_distances=False, return_indices=True)
    filled = height[tuple(nearest)]
    smooth = gaussian_filter(filled, sigma=1.0, mode='nearest')
    gx = -12.0 + (np.arange(N) + 0.5) * BIN
    gy = -12.0 + (np.arange(N) + 0.5) * BIN
    X, Y = np.meshgrid(gx, gy)
    mask = (np.hypot(X, Y) > 10.5) | (distance_transform_edt(~valid) > 2.0)
    Z = np.where(mask, np.nan, smooth)
    atom_indices = np.array(list(atoms.values()), dtype=int)
    return dx, dy, X, Y, Z, atom_indices

paths = [
    ('No-impact control', ROOT / 'structures/nohit_post_cool_0p05ps.data'),
    ('Eight 62 eV Ag impacts', ROOT / 'structures/impact_post_cool_0p05ps.data'),
]
fig = plt.figure(figsize=(13, 6.4), constrained_layout=True)
mesh = None
for col, (title, path) in enumerate(paths, 1):
    ids, types, xyz = read_data(path)
    dx, dy, X, Y, Z, atom_indices = surface_envelope(ids, xyz)
    ax = fig.add_subplot(1, 2, col, projection='3d')
    mesh = ax.plot_surface(X, Y, Z, cmap='viridis', vmin=30, vmax=52,
                           linewidth=0, antialiased=True, alpha=0.88, shade=True)
    ax.scatter(dx[atom_indices], dy[atom_indices], xyz[atom_indices, 2],
               c=xyz[atom_indices, 2], cmap='viridis', vmin=30, vmax=52,
               s=14, depthshade=False, edgecolors='black', linewidths=0.18)
    if col == 2:
        proj = ids > 4305
        px, py = local_xy(xyz[proj])
        ax.scatter(px, py, xyz[proj, 2], marker='D', s=34, color='#e23b32',
                   edgecolors='black', linewidths=0.4, label='Deposited Ag atoms')
        ax.legend(loc='upper left', fontsize=8)
    theta = np.linspace(0, 2*np.pi, 240)
    ax.plot(8*np.cos(theta), 8*np.sin(theta), np.full_like(theta, 30.0),
            color='white', linestyle='--', linewidth=1.2, alpha=0.95)
    ax.scatter([0], [0], [30], marker='x', color='red', s=55, linewidths=2)
    ax.set_xlim(-12, 12); ax.set_ylim(-12, 12); ax.set_zlim(30, 52)
    ax.set_box_aspect((24, 24, 22))
    ax.view_init(elev=33, azim=-55)
    ax.set_xlabel('Local x (Å)', labelpad=3)
    ax.set_ylabel('Local y (Å)', labelpad=3)
    ax.set_zlabel('Surface height z (Å)', labelpad=4)
    ax.set_title(title, pad=8)
    ax.grid(False)
fig.suptitle('3D surface morphology of the Ag–Ti₃SiC₂ model after 8 impacts (1.85 ps)\nSmoothed envelope from the uppermost substrate atoms; dashed ring: r = 8 Å', fontsize=13)
fig.colorbar(mesh, ax=fig.axes, shrink=0.66, pad=0.02, label='Height z (Å)')
out = ROOT / 'crater_3d_surface_8hits_1p85ps.png'
fig.savefig(out, dpi=200, bbox_inches='tight')
print(f'{out} {out.stat().st_size} bytes')
