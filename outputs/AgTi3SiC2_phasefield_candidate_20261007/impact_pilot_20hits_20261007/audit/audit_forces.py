from pathlib import Path
import numpy as np
from audit_crater import read_data, read_last_dump, min_image_xy
root=Path('.')
start=read_data('structures/cool_04_quench_0p5ps.data')
step,cols,frame=read_last_dump('dumps/cool_05_quench_0p5ps.dump')
ids=frame[:,cols.index('id')].astype(int)
xyz=frame[:,[cols.index('x'),cols.index('y'),cols.index('z')]]
forces=frame[:,[cols.index('fx'),cols.index('fy'),cols.index('fz')]]
fn=np.linalg.norm(forces,axis=1)
initial={int(i):(int(t),p) for i,t,p in zip(start['ids'],start['typ'],start['xyz'])}
fixed=set()
for line in Path('inputs/cool_05_quench_0p5ps.lmp').read_text().splitlines():
    if line.startswith('group fixed_part_'):
        fixed.update(map(int,line.split()[3:]))
r=min_image_xy(start)
freezone=set(int(i) for i,rr,p in zip(start['ids'],r,start['xyz']) if rr<8.0 and p[2]>=12 and int(i) not in fixed)
for name,sel in [('all',np.ones(len(ids),bool)),('freezone',np.array([int(i) in freezone for i in ids])),('outside_freezone',np.array([int(i) not in freezone and int(i) not in fixed for i in ids])),('Ag',frame[:,cols.index('type')]==4),('TSC',frame[:,cols.index('type')]!=4)]:
    ix=np.flatnonzero(sel)
    if len(ix):
        k=ix[np.argmax(fn[ix])]
        print(f'{name}: n={len(ix)} max_total_force={fn[k]:.6g} eV/A id={ids[k]} type={int(frame[k,cols.index("type")])} xyz={xyz[k]}')
