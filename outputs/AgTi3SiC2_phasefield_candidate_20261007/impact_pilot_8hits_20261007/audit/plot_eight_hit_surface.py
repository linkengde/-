from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
LX,LY,XY=37.28597486102,32.290601434511,-18.64298743051
CX,CY=9.9573898076,24.2578814220
CV=CY/LY; CU=((CX-XY*CV)/LX)%1

def read(path):
    lines=Path(path).read_text().splitlines(); i=next(j for j,s in enumerate(lines) if s.strip().startswith('Atoms'))+1; rows=[]
    for line in lines[i:]:
        p=line.split('#',1)[0].split()
        if len(p)<5:
            if rows: break
            continue
        try: rows.append((int(p[0]),int(p[1]),float(p[2]),float(p[3]),float(p[4])))
        except ValueError:
            if rows: break
    a=np.asarray(rows);a=a[np.argsort(a[:,0])]
    return a[:,0].astype(int),a[:,1].astype(int),a[:,2:5]

def local_xy(x):
    v=(x[:,1]/LY)%1;u=((x[:,0]-XY*(x[:,1]/LY))/LX)%1
    du=(u-CU+.5)%1-.5;dv=(v-CV+.5)%1-.5
    dx=du*LX+dv*XY;dy=dv*LY
    return dx,dy

colors={1:'#7b4ab5',2:'#4f9ca9',3:'#333333',4:'#e0b34d'}
labels={1:'Ti',2:'Si',3:'C',4:'Ag'}
paths={'No-hit control (1.8 ps)':'structures/nohit_8shot_final_cool_0p20ps.data','Eight 62 eV Ag hits (1.8 ps)':'structures/impact_8hits_final_cool_0p20ps.data'}
fig=plt.figure(figsize=(15,5.2),constrained_layout=True)
gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.05])
for col,(title,path) in enumerate(paths.items()):
    ax=fig.add_subplot(gs[0,col]);ids,t,x=read(path);dx,dy=local_xy(x)
    sel=(np.abs(dy)<2.5)&(np.abs(dx)<10)&(x[:,2]>30)
    for typ in (1,2,3,4):
        m=sel&(t==typ)&(ids<=4305)
        if m.any():ax.scatter(dx[m],x[m,2],s=18,c=colors[typ],label=labels[typ],alpha=.84,linewidths=0)
    m=sel&(ids>4305)
    if m.any():ax.scatter(dx[m],x[m,2],s=34,c='#d62728',marker='D',label='Ag impactors',edgecolors='black',linewidths=.3)
    ax.axvline(0,color='#777777',lw=.7,ls=':')
    ax.set_xlim(-10,10);ax.set_ylim(30,54)
    ax.set_xlabel('Lateral distance from impact center (Å)')
    if col==0:ax.set_ylabel('Height z (Å)')
    ax.set_title(title)
    ax.grid(alpha=.15)
    if col==0:ax.legend(fontsize=8,loc='lower left',ncol=2,frameon=True)

ax=fig.add_subplot(gs[0,2]);ids,t,x=read(paths['Eight 62 eV Ag hits (1.8 ps)']);dx,dy=local_xy(x);m=(x[:,2]>35)&(ids<=4305)
sc=ax.scatter(dx[m],dy[m],c=x[m,2],s=21,cmap='viridis',vmin=38,vmax=51,linewidths=0)
for r in (2,4,8):ax.add_patch(plt.Circle((0,0),r,fill=False,color='white' if r<5 else 'black',lw=.9,ls='--',alpha=.85))
ax.scatter([0],[0],s=85,marker='x',color='red',linewidths=2,label='impact center')
ax.set_aspect('equal');ax.set_xlim(-10,10);ax.set_ylim(-10,10);ax.set_xlabel('In-plane distance x (Å)');ax.set_ylabel('In-plane distance y (Å)');ax.set_title('Eight-hit top view\n(color = height z; substrate atoms z > 35 Å)');ax.legend(fontsize=8,loc='upper right')
fig.colorbar(sc,ax=ax,label='Height z (Å)',shrink=.78)
fig.suptitle('Ag–Ti₃SiC₂ small-cell 8-hit pilot: surface crater and topography',fontsize=14)
out=Path('audit/surface_morphology_8hits_1p80ps.png');fig.savefig(out,dpi=180,bbox_inches='tight');print(out, out.stat().st_size)
