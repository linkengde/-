from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_crater import read_data, CENTER

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / 'checkpoints/after_hit_446.data'
OUT = ROOT / 'images/crater_after446_100ps_3d.png'


def surface_grid(atoms, nbin=96):
    u = np.mod(atoms['frac'][:, 0], 1.0)
    v = np.mod(atoms['frac'][:, 1], 1.0)
    z = atoms['xyz'][:, 2]
    grid = np.full((nbin, nbin), np.nan)
    for i, j, zz in zip((u * nbin).astype(int), (v * nbin).astype(int), z):
        if zz < 0 or zz > 55:
            continue
        if not np.isfinite(grid[j, i]) or zz > grid[j, i]:
            grid[j, i] = zz
    valid = np.isfinite(grid)
    # Normalized, periodic smoothing fills empty bins without biasing the top down.
    num = gaussian_filter(np.where(valid, grid, 0.0), 4.0, mode='wrap')
    den = gaussian_filter(valid.astype(float), 4.0, mode='wrap')
    smooth = np.full_like(grid, np.nan)
    good = den > 0.08
    smooth[good] = num[good] / den[good]
    return smooth


def nearest_image(xy, atoms):
    a = np.array([atoms['lx'], 0.0])
    b = np.array([atoms['tilt'][0], atoms['ly']])
    delta = xy - CENTER
    candidates = np.stack([delta + i * a + j * b
                           for i in (-2, -1, 0, 1, 2)
                           for j in (-2, -1, 0, 1, 2)])
    idx = np.argmin(np.sum(candidates * candidates, axis=2), axis=0)
    return candidates[idx, np.arange(len(xy))]


def make_local_mesh(atoms, grid, nbin=96):
    us = (np.arange(nbin) + 0.5) / nbin
    vs = (np.arange(nbin) + 0.5) / nbin
    vv, uu = np.meshgrid(vs, us, indexing='ij')
    x = atoms['bounds']['x'][0] + uu * atoms['lx'] + vv * atoms['tilt'][0]
    y = atoms['bounds']['y'][0] + vv * atoms['ly']
    xy = np.column_stack((x.ravel(), y.ravel()))
    local = nearest_image(xy, atoms).reshape(nbin, nbin, 2)
    radius = np.linalg.norm(local, axis=2)
    z = np.ma.masked_where((radius > 18.0) | ~np.isfinite(grid), grid)
    return local[:, :, 0], local[:, :, 1], z


def main():
    final = read_data(FINAL)
    grid = surface_grid(final)
    X, Y, Z = make_local_mesh(final, grid)

    # Show the actual atomistic surface shell. Airborne/sputtered atoms above z=55 Å
    # are left out only to keep the crater and the two phases legible.
    xyz = final['xyz']
    frac = final['frac']
    u = np.mod(frac[:, 0], 1.0)
    v = np.mod(frac[:, 1], 1.0)
    ii = np.minimum((u * grid.shape[1]).astype(int), grid.shape[1] - 1)
    jj = np.minimum((v * grid.shape[0]).astype(int), grid.shape[0] - 1)
    local_top = grid[jj, ii]
    pos = nearest_image(xyz[:, :2], final)
    radius = np.linalg.norm(pos, axis=1)
    z = xyz[:, 2]
    shell = (radius < 18.0) & (z > 0.0) & (z < 55.0) & np.isfinite(local_top) & (z >= local_top - 5.0)

    plt.rcParams.update({'font.size': 11, 'axes.titleweight': 'bold'})
    fig = plt.figure(figsize=(12.5, 9.0), dpi=190)
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='inferno', vmin=18, vmax=43,
                           linewidth=0, antialiased=True, alpha=0.86,
                           rcount=90, ccount=90, shade=True)
    # Split the actual near-surface atoms by phase: purple Ti3SiC2, gold Ag.
    colors = {1: '#7851a9', 2: '#a778cf', 3: '#d0b6e8', 4: '#e4b647'}
    sizes = {1: 7, 2: 7, 3: 7, 4: 8}
    labels = {1: 'Ti', 2: 'Si', 3: 'C', 4: 'Ag'}
    for typ in (1, 2, 3, 4):
        m = shell & (final['typ'] == typ)
        ax.scatter(pos[m, 0], pos[m, 1], z[m], s=sizes[typ], c=colors[typ],
                   alpha=0.96, edgecolors='none', depthshade=False, label=labels[typ])
    fcenter = np.linalg.solve(final['H'], CENTER - np.array([final['bounds']['x'][0], final['bounds']['y'][0]])) % 1.0
    jc = min(int(fcenter[1] * grid.shape[0]), grid.shape[0] - 1)
    ic = min(int(fcenter[0] * grid.shape[1]), grid.shape[1] - 1)
    center_z = grid[jc, ic] if np.isfinite(grid[jc, ic]) else 24.0
    ax.scatter([0], [0], [center_z], s=36, marker='x', color='cyan', linewidths=2.2,
               label='Impact centre')
    ax.set_xlim(-18, 18)
    ax.set_ylim(-18, 18)
    ax.set_zlim(12, 48)
    ax.set_xlabel('x from impact centre (Å)', labelpad=10)
    ax.set_ylabel('y from impact centre (Å)', labelpad=10)
    ax.set_zlabel('z (Å)', labelpad=7)
    ax.set_title('Ag–Ti₃SiC₂ after 100 ps cumulative Ag bombardment\n3D top surface with near-surface atoms', pad=18, fontsize=15)
    ax.view_init(elev=31, azim=-48)
    ax.set_box_aspect((36, 36, 28))
    ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.96), ncol=2,
              framealpha=0.9, fontsize=9)
    cb = fig.colorbar(surf, ax=ax, shrink=0.62, pad=0.08, aspect=24)
    cb.set_label('smoothed top-surface height z (Å)')
    fig.text(0.5, 0.035,
             'Surface reconstructed from the final atom coordinates; Ag/Ti₃SiC₂ atoms shown within 5 Å of it. '
             '128 atoms above z=55 Å omitted from view. Height mesh is smoothed for visibility.',
             ha='center', va='center', fontsize=8.7)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(OUT)
    print(f'{OUT.stat().st_size} bytes')


if __name__ == '__main__':
    main()
