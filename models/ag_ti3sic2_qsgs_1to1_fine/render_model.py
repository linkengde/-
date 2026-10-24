#!/usr/bin/env python3
"""Scientific morphology previews, explicitly distinguished from atomic files."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from skimage.measure import marching_cubes
import plotly.graph_objects as go

ROOT = Path(__file__).resolve().parent
data = np.load(ROOT/'qsgs_morphology.npz')
mask = data['skeleton']
size = data['size_A']
spacing = size / np.array(mask.shape)
cmap = ListedColormap(['#bcbfc4', '#205dcc'])
fig, axes = plt.subplots(1, 3, figsize=(13, 5))
for ax, array, extent, label, xlabel, ylabel in [
    (axes[0], mask[:, :, 25].T, [0, size[0], 0, size[1]], f'XY at z = {25.5*spacing[2]:.1f} A', 'x (A)', 'y (A)'),
    (axes[1], mask[:, 50, :].T, [0, size[0], 0, size[2]], f'XZ at y = {50.5*spacing[1]:.1f} A', 'x (A)', 'z (A)'),
    (axes[2], mask[50, :, :].T, [0, size[1], 0, size[2]], f'YZ at x = {50.5*spacing[0]:.1f} A', 'y (A)', 'z (A)')]:
    ax.imshow(array, origin='lower', interpolation='nearest', cmap=cmap, vmin=0, vmax=1, extent=extent)
    ax.set(title=label, xlabel=xlabel, ylabel=ylabel)
fig.suptitle(f'QSGS target morphology | Blue: Ti3SiC2; Grey: Ag domain\n{size[0]:.2f} x {size[1]:.2f} x {size[2]:.2f} A | Atomic equal-mass candidate; unrelaxed', fontsize=12)
fig.tight_layout()
fig.savefig(ROOT/'morphology_sections.png', dpi=180)
plt.close(fig)

# Expose the actual solid domains at the exterior by zero padding.
vertices, faces, _, _ = marching_cubes(np.pad(mask.astype(np.float32), 1), .5, spacing=tuple(spacing), step_size=2)
vertices -= spacing/2
vertices = np.clip(vertices, [0, 0, 0], size)
mesh = go.Mesh3d(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2], i=faces[:, 0], j=faces[:, 1], k=faces[:, 2], color='#205dcc', flatshading=True, name='Ti3SiC2 QSGS domain', opacity=1)
silver_vertices, silver_faces, _, _ = marching_cubes(np.pad((~mask).astype(np.float32), 1), .5, spacing=tuple(spacing), step_size=2)
silver_vertices = np.clip(silver_vertices-spacing/2, [0, 0, 0], size)
ag = go.Mesh3d(x=silver_vertices[:, 0], y=silver_vertices[:, 1], z=silver_vertices[:, 2], i=silver_faces[:, 0], j=silver_faces[:, 1], k=silver_faces[:, 2], color='#c3c6cb', flatshading=True, name='Ag complementary domain', opacity=.18)
fig = go.Figure([mesh, ag])
fig.update_layout(title='QSGS morphology preview (coarsened voxel surface, not atomic relaxation)<br>Blue: Ti3SiC2; toggle Ag in legend | Ag:Ti3SiC2 atomic mass = 1:1', scene=dict(xaxis_title='x (A)', yaxis_title='y (A)', zaxis_title='z (A)', aspectmode='data'), showlegend=True, margin=dict(l=0,r=0,b=0,t=70))
fig.write_html(ROOT/'morphology_3d.html', include_plotlyjs=True, full_html=True)

# Exact final atomic occupancy slices, including interface clearance voids.
from collections import defaultdict
with (ROOT/'ag_ti3sic2_qsgs.extxyz').open() as handle:
    n = int(handle.readline());handle.readline()
    # Avoid a second dependency on ASE in this plotting-only environment.
    pos = []; species = []
    for line in handle:
        f=line.split();species.append(f[0]);pos.append([float(x) for x in f[1:4]])
pos=np.array(pos);species=np.array(species)
fig,ax=plt.subplots(figsize=(7,7))
mid_z = size[2]/2
sel=np.abs(pos[:,2]-mid_z)<1.5
order=np.argsort(species[sel]=='Ag')
xy=pos[sel,:2][order];sp=species[sel][order]
ax.scatter(xy[:,0],xy[:,1],c=np.where(sp=='Ag','#b9bbc1','#205dcc'),s=1.8,linewidths=0)
ax.set(xlim=(0,size[0]),ylim=(0,size[1]),aspect='equal',xlabel='x (A)',ylabel='y (A)',title=f'Final atomic slice: z = {mid_z:.2f} +/- 1.5 A\nBlue: Ti3SiC2 atoms; grey: Ag; white: unoccupied')
fig.tight_layout();fig.savefig(ROOT/'atomic_xy_slice.png',dpi=180)
print('Saved morphology_sections.png, morphology_3d.html and atomic_xy_slice.png')
