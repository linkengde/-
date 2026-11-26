import numpy as np, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from pathlib import Path
R=Path(__file__).parent
font=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')

def load_ag_dump(path):
    with open(path) as f:
        while 1:
            l=f.readline()
            if not l: raise RuntimeError(path)
            if l.startswith('ITEM: TIMESTEP'):
                step=int(f.readline()); f.readline(); n=int(f.readline()); f.readline()
                for _ in range(3): f.readline()
                h=f.readline().split()[2:]; c={x:i for i,x in enumerate(h)}
                a=np.empty((n,len(h)))
                for i in range(n): a[i]=np.fromstring(f.readline(),sep=' ')
                a=a[a[:,c['type']]==4]
                return {k:a[:,c[k]] for k in ('x','y','z')}

def radial_surface(d):
    x=d['x']-100.; y=d['y']-100.; z=d['z']; k=(np.abs(x)<50)&(np.abs(y)<50)&(z>=150)&(z<=230)
    tiles={}
    for xx,yy,zz in zip(x[k],y[k],z[k]):
        key=(int(np.floor((xx+50)/4)),int(np.floor((yy+50)/4)))
        tiles[key]=max(tiles.get(key,-1e9),float(zz))
    r=[]; h=[]
    for (ix,iy),zz in tiles.items(): r.append(np.hypot((ix+.5)*4-50,(iy+.5)*4-50)); h.append(zz)
    r=np.array(r); h=np.array(h)
    centers=np.arange(5,51,10); med=[]; q25=[]; q75=[]
    for center in centers:
        m=(r>=center-5)&(r<center+5)
        med.append(float(np.median(h[m])) if np.any(m) else np.nan)
        q25.append(float(np.percentile(h[m],25)) if np.any(m) else np.nan)
        q75.append(float(np.percentile(h[m],75)) if np.any(m) else np.nan)
    return centers,np.array(med),np.array(q25),np.array(q75)

pre=load_ag_dump(R/'seed41_local_pulse_1479eV_start.dump')
post=load_ag_dump(R/'seed41_local_pulse_final_5p8ps_all.dump')
r0,m0,l0,h0=radial_surface(pre); r1,m1,l1,h1=radial_surface(post)
rows=list(csv.DictReader(open(R/'seed41_local_pulse_frame_metrics.csv')))
t=np.array([float(x['time_ps']) for x in rows]); n220=np.array([int(x['n_z_gt_220']) for x in rows]); n250=np.array([int(x['n_z_gt_250']) for x in rows]); n300=np.array([int(x['n_z_gt_300']) for x in rows]); local=np.array([int(x['n_local_r25_z_gt_220']) for x in rows])
fig,ax=plt.subplots(1,2,figsize=(11,4.7),constrained_layout=True)
ax[0].fill_between(r0,l0,h0,color='#8ecae6',alpha=.25); ax[0].plot(r0,m0,'o-',color='#4d95b5',label='脉冲前 t=1.8 ps')
ax[0].fill_between(r1,l1,h1,color='#263f4e',alpha=.14); ax[0].plot(r1,m1,'s-',color='#263f4e',label='脉冲后 t=5.8 ps')
ax[0].axvspan(0,25,color='#ff8c24',alpha=.08,label='脉冲半径 25 Å')
ax[0].set_xlabel('距脉冲中心的径向距离 r (Å)',fontproperties=font); ax[0].set_ylabel('Ag 表面高度中位数 z (Å)',fontproperties=font)
ax[0].set_title('表面径向轮廓',fontproperties=font); ax[0].grid(alpha=.25); ax[0].legend(prop=font,fontsize=9)
ax[1].plot(t,n220,'o-',ms=3,color='#4d95b5',label='全局 Ag：z > 220 Å')
ax[1].plot(t,n250,'s-',ms=3,color='#e76f51',label='全局 Ag：z > 250 Å')
ax[1].plot(t,local,'^-',ms=3,color='#f4a261',label='中心 r < 25 Å：z > 220 Å')
ax[1].plot(t,n300,'--',color='#6a4c93',label='全局 Ag：z > 300 Å')
ax[1].set_xlabel('时间 t (ps)',fontproperties=font); ax[1].set_ylabel('Ag 原子数',fontproperties=font)
ax[1].set_title('进入上方真空区的 Ag 计数',fontproperties=font); ax[1].grid(alpha=.25); ax[1].legend(prop=font,fontsize=9)
fig.suptitle('seed41：局域 1.479 keV 脉冲后的表面变化与 Ag 喷出统计',fontproperties=font,fontsize=14)
fig.text(.5,-.035,'中心 r<10 Å 表面中位高度下降约 4.14 Å；中心相对 r=20–40 Å 周围额外降低约 5.03 Å，中心 Ag 数减少 106 个，支持浅凹坑趋势但仍需谨慎。末态：209 个 Ag 在 z>220 Å、43 个在 z>250 Å、6 个在 z>300 Å；均未离开 500 Å 计算盒。',ha='center',fontproperties=font,fontsize=9)
out=R/'seed41_crater_ejection_trends_300dpi.png'; fig.savefig(out,dpi=300,bbox_inches='tight',facecolor='white'); print(out,out.stat().st_size)
