"""使用已审核的原面CSV和原生诊断绘科研图；不运行求解或改变工作台默认结果。"""
from pathlib import Path
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def main():
    start=time.monotonic();root=Path(__file__).resolve().parents[1];base=root/'outputs/v197_validation'
    r=json.loads((base/'rocket_audit01/report.json').read_text());case=base/'rocket01'
    s=json.loads((case/'surface.json').read_text());v=np.array(s['vertices']);tri=v[np.array(s['facesOneBased'])-1]
    wall=np.loadtxt(base/'rocket_audit01/wall_diagnostic.csv',delimiter=',',skiprows=1)
    names=(case/'lmrsim/diagnostics/nk-diagnostics').read_text().splitlines()[0].split()
    raw=np.loadtxt(case/'lmrsim/diagnostics/nk-diagnostics',skiprows=1,ndmin=2)
    history={name:raw[:,i] for i,name in enumerate(names)}
    fig=plt.figure(figsize=(13,8.5));status='Development screening passed; NOT converged' if r['candidateScreeningPassed'] else 'Candidate rejected; diagnostic values only'
    fig.suptitle('Original 19,752-face rocket | M5, alpha 20 deg | '+status,fontsize=12)
    lo,hi=v.min(0),v.max(0)
    for index,column,label,limits in [(1,7,'Fluid wall-neighbor T (K), NOT wall T',(200,2200)),
                                      (2,6,'Cell-average pressure (Pa)',(0,40000))]:
        ax=fig.add_subplot(2,2,index,projection='3d');norm=matplotlib.colors.Normalize(*limits)
        colors=plt.cm.cividis(norm(wall[:,column]));colors[wall[:,8]==0]=[.85,.16,.2,1]
        ax.add_collection3d(Poly3DCollection(tri,facecolors=colors,linewidth=0))
        ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]))
        ax.set_box_aspect(hi-lo);ax.view_init(20,-70);ax.set_axis_off()
        ax.set_title('Saved Newton step %d | fixed color limits'%r['savedSteps'][-1],fontsize=10,y=.94)
        ax.text2D(.02,.05,'Full original area: 1.9646292566 m2\nHeat flux and solved solid wall T: unknown',transform=ax.transAxes,fontsize=8)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='cividis'),ax=ax,shrink=.65,pad=0,label=label,extend='both')
    ax=fig.add_subplot(223)
    for col,label in [('mass-rel','mass'),('x-mom-rel','x momentum'),('y-mom-rel','y momentum'),('z-mom-rel','z momentum'),('total-energy-rel','energy')]:
        ax.semilogy(history['step'],history[col],label=label)
    ax.axhline(1e-8,color='gray',ls='--',label='Required 1e-8')
    ax.set(xlabel='Newton update',ylabel='Native component / initial reference',title='Rocket residuals: numerical gate remains independent')
    ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=2)
    ax=fig.add_subplot(224)
    control=np.loadtxt(base/'audit01/wall_development.csv',delimiter=',',skiprows=1)
    ref=json.loads((base/'audit01/same_eos_reference.json').read_text())
    ax.scatter(-control[:,1],control[:,5]/1000,s=12,alpha=.5,label='14,400-tet control, step 120')
    for i,span in [(1,(.35,.5)),(2,(.85,.95))]:
        ax.hlines(ref['states'][i]['pressurePa']/1000,*span,color='black',ls='--',label='Same-EOS inviscid reference' if i==1 else None)
    ax.set(xlabel='Control axial distance -X (m)',ylabel='Wall-neighbor pressure (kPa)',title='Small control: developing, not accuracy-qualified')
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout(rect=(0,0,1,.96));path=base/'nk_development.png';fig.savefig(path,dpi=180);plt.close(fig)
    print(json.dumps(dict(passed=True,path=str(path),elapsedSeconds=time.monotonic()-start)))


if __name__=='__main__':main()
