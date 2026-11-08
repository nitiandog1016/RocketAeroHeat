"""用真实原面和已有保存场绘图；范围示意不冒充压力/热流云图。"""
from pathlib import Path
import json,time,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from prepareOpenFoamControl import ROOT,save,sha

def main():
    start=time.monotonic();base=ROOT/'outputs/v207_validation';out=base/'figures01';out.mkdir(exist_ok=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    s=json.loads((ROOT/'outputs/v194_validation/rocket_design03/surface.json').read_text())
    c=json.loads((base/'region01/contract.json').read_text());p=np.array(s['vertices']);f=np.array(s['facesOneBased'])-1
    ids=np.array(c['wallOriginalFaceIDs'])-1;selected=np.zeros(len(f),bool);selected[ids]=True
    data=np.load(base/'region01/boundary.npz');v=data['points'];tri=data['faces'];patch=data['patch']
    # 模块一：整箭保持完整分母；局部仅以单色高亮准备范围。
    fig=plt.figure(figsize=(13,5.4));a=fig.add_subplot(121,projection='3d');b=fig.add_subplot(122,projection='3d')
    a.add_collection3d(Poly3DCollection(p[f[~selected]],facecolor='#c4c9ce',linewidth=0))
    a.add_collection3d(Poly3DCollection(p[f[selected]],facecolor='#267d9c',linewidth=0))
    a.set(xlim=(p[:,0].min(),p[:,0].max()),ylim=(-.5,.5),zlim=(-.5,.5),xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)')
    a.set_box_aspect((np.ptp(p[:,0]),1,1));a.view_init(elev=22,azim=-60)
    a.set_title('Original rocket: prepared wall region',pad=10)
    a.text2D(.02,.02,'7,838 / 19,752 original faces\n0.719860 / 1.964629 m$^2$ (36.64%)',transform=a.transAxes)
    b.add_collection3d(Poly3DCollection(v[tri[patch==0]],facecolor='#267d9c',linewidth=0))
    # 画外表面及环形出口，原物理壁没有人工实体端盖。
    for tag,color in [(1,'#7ca59b'),(2,'#819ab3'),(3,'#bd9355')]:
        b.add_collection3d(Poly3DCollection(v[tri[patch==tag]],facecolor=color,alpha=.12,edgecolor=color,linewidth=.15))
    r=c['outletRadiusM'];b.set(xlim=(c['cutXM'],c['inletXM']),ylim=(-r,r),zlim=(-r,r),xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)')
    b.set_box_aspect((c['inletXM']-c['cutXM'],2*r,2*r));b.view_init(elev=20,azim=-55)
    b.set_title('Closed fluid boundary; volume mesh pending',pad=10)
    b.text2D(.02,.02,'Annular fluid outlet at X = -0.89446 m\nOuter extent is a design, not a verified shock envelope',transform=b.transAxes,fontsize=9)
    fig.suptitle('v207 | Geometry preparation only — no regional CFD / no new qualified heat area',fontsize=12)
    fig.subplots_adjust(left=.02,right=.98,top=.88,bottom=.08,wspace=.15)
    fig.savefig(out/'region_geometry.png',dpi=160);plt.close(fig)
    # 模块二：实际四帧热流与独立稳态BVP比较，不能据此宣称高Mach/三维收敛。
    first=json.loads((base/'thermal_audit01/report.json').read_text())
    last=json.loads((base/'thermal_audit03/report.json').read_text())
    oracle=json.loads((base/'period_repair01/oracle.json').read_text())
    fig,(a,b)=plt.subplots(1,2,figsize=(11,4.5))
    colors={'cold':'#267d9c','hot':'#b6792c'}
    for name in ['cold','hot']:
        for report,style,label in [(first,'--','upwind (failed)'),(last,'-','vanLeer (passed)')]:
            a.plot([float(x['time'])*1e9 for x in report['records']],
                   [100*x['heatFluxRelativeError'][name] for x in report['records']],style+'o',color=colors[name],label=name+' / '+label)
    a.axhline(2,color='#555555',lw=1,label='Original 2% gate');a.set(xlabel='Physical time (ns)',ylabel='Wall heat-flux relative error (%)')
    a.legend(fontsize=8);a.grid(alpha=.2)
    x=np.arange(2);b.bar(x-.18,[oracle['wallHeatIntoWallWm2'][n]/1000 for n in ['cold','hot']],.36,color='#a6adb4',label='Independent BVP')
    b.bar(x+.18,[last['records'][-1]['heatFluxIntoWallWm2'][n]/1000 for n in ['cold','hot']],.36,color='#267d9c',label='Native saved field')
    b.axhline(0,color='#555555',lw=1);b.set(xticks=x,xticklabels=['Cold wall: 300 K','Hot wall: 600 K'],ylabel='Heat into wall (kW/m$^2$)');b.legend(fontsize=9)
    fig.suptitle('v207 | Variable-property Couette interface check — 32 cells, 8 steps',fontsize=12)
    fig.text(.5,.01,'Not an experiment, a converged high-Mach solution, or new rocket heat coverage.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,.94]);fig.savefig(out/'thermal_interface.png',dpi=160);plt.close(fig)
    shutil.copyfile(__file__,out/'executed_entry.py')
    save(out/'report.json',dict(passed=True,figures={p.name:sha(p) for p in out.glob('*.png')},
        elapsedSeconds=time.monotonic()-start,physicalQualification=False))
    print(str(out))

if __name__=='__main__':main()
