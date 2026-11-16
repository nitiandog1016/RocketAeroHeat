"""模块职责：在完整原箭体上展示区域责任及真实接口特征，不把旧场画成新热解。"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from runSU2RocketPhaseV222 import ROOT,read,save,sha
import prepareViscousRegionV244 as run
BASE=run.BASE;DOC=run.DOC;CLOSED=False


def render():
    out=BASE/'visual01';out.mkdir(exist_ok=False)
    audit=read(BASE/'audit02/report.json');assert audit['passed']
    g=np.load(str(ROOT/'outputs/v230_validation/audit02/geometry.npz'))
    print('geometry keys',g.files)
    xyz=g['vertices'] if 'vertices' in g.files else g['points']
    faces=g['faces'];idx=np.load(str(run.INDEX));front=idx['upstreamOriginalFaces']-1
    mask=np.zeros(len(faces),bool);mask[front]=True
    table=np.loadtxt(str(BASE/'region02/interface.csv'),delimiter=',',skiprows=1)
    fig=plt.figure(figsize=(14,9),facecolor='white');grid=fig.add_gridspec(2,2,height_ratios=[.8,1.2],hspace=.50,wspace=.32)
    ax=fig.add_subplot(grid[0,:]);tri=xyz[faces][:,:,[0,2]]
    colors=np.tile(np.array([.76,.79,.81,1]),(len(faces),1));colors[mask]=[.13,.43,.60,1]
    ax.add_collection(PolyCollection(tri,facecolors=colors,edgecolors='none'));ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim')
    ax.set_xlabel('Original X [m]');ax.set_ylabel('Original Z [m]')
    ax.set_title('Complete original rocket / X-Z projection: 4,434 front parents (blue), 15,318 downstream parents (grey)',fontsize=11)
    ax.grid(alpha=.18);ax.margins(x=.015,y=.15)
    left=fig.add_subplot(grid[1,0]);right=fig.add_subplot(grid[1,1])
    shades={0:'#6b8796',1:'#c47f26',4:'#b64b4b',5:'#743b78'}
    radius=float(np.linalg.norm(xyz[:,1:3],axis=1).max());limit=1.4*radius
    for count,color in shades.items():
        sel=table[:,14]==count;near=sel&(abs(table[:,5])<limit)&(abs(table[:,6])<limit)
        left.scatter(table[near,5],table[near,6],s=10,c=color,label='%d incoming: %d full-interface faces'%(count,sum(sel)),alpha=.8)
        right.scatter(table[sel,17],table[sel,13],s=8,c=color,alpha=.65)
    left.set_xlim(-limit,limit);left.set_ylim(-limit,limit);left.set_aspect('equal');left.grid(alpha=.18)
    left.set_xlabel('Interface centroid Y [m]');left.set_ylabel('Interface centroid Z [m]')
    left.set_title('Near-body interface view (all 1,984 faces counted in legend)',fontsize=10)
    left.legend(loc='upper center',bbox_to_anchor=(.5,-.18),fontsize=9,ncol=2,frameon=False)
    right.set_xscale('symlog',linthresh=1e-5);right.axhline(1,c='#343a40',lw=.8,ls='--');right.axhline(0,c='#343a40',lw=.8)
    right.axhline(-1,c='#343a40',lw=.8,ls='--');right.grid(alpha=.18)
    right.set_xlabel('Minimum saved full-domain wall distance [m] / symlog')
    right.set_ylabel('Outward normal Mach from mean conserved state')
    right.set_title('Actual saved v236/300 research state / no new CFD',fontsize=10)
    fig.suptitle('RocketAeroHeat | Viscous region boundary evidence\n827 interfaces require downstream information; fixed-cut production rejected',fontsize=16,y=.98)
    fig.text(.08,.025,'Original 19,752 faces preserved. 46 interfaces touch the wall by topology; only 14 have an exactly zero stored distance.\nFace-level characteristic classification is not native dual-face flux, a proof of separation, or a qualified heating result.',fontsize=10,color='#444')
    fig.subplots_adjust(top=.84,bottom=.18,left=.08,right=.97)
    path=out/'viscous_region_boundary.png';fig.savefig(str(path),dpi=150);plt.close(fig)
    save(out/'report.json',dict(passed=True,originalFaces=19752,frontParents=4434,
        interfaceFaces=1984,sourceResearchUpdate=300,newCFDUpdates=0,physicalQualification=False,
        imageSHA256=sha(path),dataSHA256=sha(BASE/'region02/interface.csv')))


def main():
    if CLOSED:raise SystemExit('v244 region render closed; use NEXT_DECISION')
    run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC;size=run.guard.size_bytes
    run.guard.size_bytes=lambda p:size(BASE) if Path(p)==BASE or Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
        BASE/'audit02/report.json',BASE/'region02/interface.csv',run.INDEX,
        ROOT/'outputs/v230_validation/audit02/geometry.npz']
    return run.guard.supervise('render',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v244 region render closed; use NEXT_DECISION')
    if sys.argv[1]=='task':render()
    else:raise SystemExit(0 if main() else 1)
