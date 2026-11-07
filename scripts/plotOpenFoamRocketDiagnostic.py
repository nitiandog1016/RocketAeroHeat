"""原三维箭体定量诊断图：只显示首步p/流体T，明确失败任务及尚未发展的身份。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from readOpenFoamRocketState import BASE,read_json


def main():
    start=time.monotonic();out=BASE/'plot03';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    source=BASE/'rocket_audit02';r=read_json(source/'report.json');assert r['auditPassed'] and not r['planAcceptance']
    data=np.load(source/'original_wall_state.npz');s=read_json(ROOT/'outputs/v194_validation/rocket_design03/surface.json')
    v=np.array(s['vertices']);tri=v[np.array(s['facesOneBased'])-1]
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import Normalize
    fig=plt.figure(figsize=(12,4.2),facecolor='white')
    for i,(key,label,factor) in enumerate([('pressurePa','Pressure (kPa)',.001),('fluidTemperatureK','Adjacent-fluid temperature (K)',1)]):
        ax=fig.add_subplot(1,2,i+1,projection='3d');values=data[key]*factor;norm=Normalize(float(values.min()),float(values.max()))
        cmap=plt.get_cmap('viridis');collection=Poly3DCollection(tri,facecolors=cmap(norm(values)),edgecolors='none',linewidth=0)
        ax.add_collection3d(collection);ax.set_xlim(v[:,0].min(),v[:,0].max());ax.set_ylim(v[:,1].min(),v[:,1].max());ax.set_zlim(v[:,2].min(),v[:,2].max())
        ax.set_box_aspect(np.ptp(v,axis=0));ax.view_init(elev=22,azim=-64)
        ax.set_yticks([]);ax.set_zticks([])
        ax.text2D(.5,.05,'X (m); nose at +X',transform=ax.transAxes,ha='center',fontsize=9)
        ax.set_xticks([-2.5,-1.5,-.5,0]);ax.set_title(label,pad=5)
        m=plt.cm.ScalarMappable(norm=norm,cmap=cmap);m.set_array(values);fig.colorbar(m,ax=ax,shrink=.55,pad=.04)
        ax.tick_params(labelsize=8);ax.grid(False)
    fig.suptitle('Original 19,752-face rocket | Mach 5, alpha 20 deg\nFirst-step diagnostic only: t = %.6f microseconds'% (r['savedState']['nativeTimeSeconds']*1e6),fontsize=13)
    fig.text(.5,.055,'Run stopped later on a temperature-range warning. Unconverged fluid state; wall heat flux and solid temperature are unknown.',ha='center',fontsize=9,color='#7b3f2b')
    fig.subplots_adjust(left=.01,right=.96,top=.78,bottom=.13,wspace=.13);fig.savefig(out/'rocket_first_step.png',dpi=190);plt.close(fig)
    save(out/'report.json',dict(passed=True,sourceReportSHA256=sha(source/'report.json'),sourceStateSHA256=sha(source/'original_wall_state.npz'),
        originalWallFaces=19752,scientificMeaning='first-step diagnostic; not developed pressure/heat solution',
        mainUIChanged=False,CFDRuns=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start))
    print(json.dumps(dict(passed=True,figure=str(out/'rocket_first_step.png'),seconds=time.monotonic()-start)))


if __name__=='__main__':main()
