# -*- coding: utf-8 -*-
"""比较同几何已收敛粗/细Euler壁邻解；按面积与同一轴向区间，不混用节点编号。"""
from pathlib import Path
import argparse
import json
import shutil
import time
import numpy as np
from readHybridEulerNative import read_wall_frame
from readEilmerNative import read_grid,read_field,sha


def compare(base,out):
    base,out=map(Path,(base,out));out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    shutil.copyfile(__file__,out/'executed_comparison.py')
    paths=[base/'coarse_audit02/report.json',base/'fine_audit02/report.json'];reports=[json.loads(p.read_text()) for p in paths]
    assert all(r['passed'] for r in reports)
    cases=[base/'coarse01',base/'fine01'];preps=[base/'prepare01/report.json',base/'fine_prepare01/report.json']
    runs=[base/'coarse_run01/report.json',base/'fine_run01/report.json'];frames=[];curves=[];boundaries=[]
    refs=json.loads((base.parent/'v191_validation/audit01/reference.json').read_text())
    for name,case,prep,run,report in zip(['coarse','fine'],cases,preps,runs,reports):
        snapshot=len(report['savedSteps'])-1
        frame,info=read_wall_frame(case,snapshot,prep,run);frames.append(frame)
        (out/(name+'_final_wall.json')).write_text(json.dumps(frame,indent=2)+'\n')
        # 相同粗轴向区间上积分；细壁邻点距壁更近，解释中保留采样误差来源。
        s=-np.array(frame['faceCentersM'])[:,0];area=np.array(frame['faceAreasM2'])
        group=np.floor(s*60).astype(int);assert group.min()==0 and group.max()==59
        curve=[]
        for i in range(60):
            mask=group==i;w=area[mask]
            curve.append([float(w.sum())]+[float(np.average(np.asarray(frame['fields'][key])[mask],weights=w))
                         for key in ['pressurePa','temperatureK']])
        curves.append(np.array(curve))
        grid=read_grid(case/'lmrsim/grid/grid-0000');raw=read_field(case,snapshot,info['cells'])
        velocity=np.column_stack([raw['vel.'+k] for k in 'xyz']);binfo={}
        for key,b in grid['boundaries'].items():
            vertices=grid['points'][grid['faces'][b['ids']]]
            vector=(np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0])+
                    np.cross(vertices[:,2]-vertices[:,0],vertices[:,3]-vertices[:,0]))/2
            vector*=b['signs'][:,None];n=vector/np.linalg.norm(vector,axis=1)[:,None];idx=b['owners']
            normal=np.sum(velocity[idx]*n,axis=1)/raw['a'][idx]
            binfo[key]=dict(minOutwardNormalMach=float(normal.min()),maxOutwardNormalMach=float(normal.max()),
                           maximumPressureRelativeToFreestreamDeviation=float(np.max(abs(raw['p'][idx]/5000-1))))
        boundaries.append(binfo)
    # 独立均匀几何真值检查，避免在两份读回面积相等时忽略共同分母错误。
    exact=.04*(.6+.4/np.cos(np.radians(10)))
    assert all(abs(curve[:,0].sum()-exact)<1e-14 for curve in curves)
    assert np.max(abs(curves[0][:,0]-curves[1][:,0]))<1e-14
    difference={key:float(np.sum(curves[1][:,0]*abs(curves[1][:,j]-curves[0][:,j]))/
                           np.sum(curves[1][:,0]*abs(curves[1][:,j]))) for j,key in [(1,'pressurePa'),(2,'temperatureK')]}
    platforms=[]
    for a,b in zip(reports[0]['platformComparisons'],reports[1]['platformComparisons']):
        assert a['region']==b['region']
        platforms.append(dict(region=a['region'],coarsePressureError=a['pressureRelativeError'],finePressureError=b['pressureRelativeError'],
            coarseTemperatureError=a['temperatureRelativeError'],fineTemperatureError=b['temperatureRelativeError'],
            fineDirectionErrorDeg=b['directionErrorDeg'],
            pressureErrorReduced=abs(b['pressureRelativeError'])<abs(a['pressureRelativeError']),
            temperatureErrorReduced=abs(b['temperatureRelativeError'])<abs(a['temperatureRelativeError'])))
    result=dict(passed=True,bothNumericallyConverged=all(r['numericalConvergencePassed'] for r in reports),
        originalCoarseSamplingFailureRetained=not reports[0]['scheduledSamplingPassed'],fineSamplingPassed=reports[1]['scheduledSamplingPassed'],
        wallAreaM2=exact,areaL1FineVersusCoarseOnCommonStrips=difference,platforms=platforms,boundaryAdjacentStates=boundaries,
        physicalQualification=False,gridIndependenceProven=False,wallEdgeQualified=False,newSolvedHeatAreaM2=0,
        note='two-grid wall-adjacent comparison; changes combine spatial discretization and sampling distance; no external-domain sensitivity proof',
        sourceSHA256={str(p):sha(p) for p in paths+preps+runs+[Path(__file__)]})
    np.savetxt(out/'common_axial_strips.csv',np.column_stack([(np.arange(60)+.5)/60,curves[0],curves[1][:,1:]]),delimiter=',',
        header='s_mid_m,area_m2,coarse_p_Pa,coarse_T_K,fine_p_Pa,fine_T_K',comments='',fmt='%.17g')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,6.6),constrained_layout=True)
    for ax,col,field,label in [(axes[0,0],1,'pressurePa','Wall-adjacent pressure [Pa]'),(axes[0,1],2,'temperatureK','Temperature [K]')]:
        for curve,name in zip(curves,['60 x 2 x 20','120 x 2 x 40']):ax.plot((np.arange(60)+.5)/60,curve[:,col],label=name,lw=1.4)
        for lo,hi,ref in [(.35,.5,refs['states'][1]),(.85,.95,refs['states'][2])]:ax.plot([lo,hi],[ref[field]]*2,'k--',lw=1.2)
        ax.set(xlabel='s = -X [m]',ylabel=label);ax.legend(fontsize=8);ax.grid(alpha=.2)
    for ax,case,label in zip(axes[1],cases,['Coarse restart: step 400 to 2000','Fine: step 0 to 2400']):
        history=np.loadtxt(case/'lmrsim/transient-residuals.dat');rho=refs['states'][0]['densityKgm3'];U=np.linalg.norm(refs['states'][0]['velocityMps'])
        ax.semilogy(history[:,0],history[:,3]/(rho*U),label='mass')
        ax.semilogy(history[:,0],history[:,7]/(rho*U**3),label='energy')
        ax.axhline(1e-8,color='k',ls='--',lw=1);ax.set(xlabel='Pseudo-time update',ylabel='Dimensionless Linf residual',title=label)
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Same-EOS Euler compression / expansion: iteration vs. spatial error\nTwo grids; no viscous heat flux or rocket physical qualification',fontsize=11)
    fig.savefig(out/'refinement_comparison.png',dpi=160);plt.close(fig)
    result['elapsedSeconds']=time.monotonic()-start
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();compare(a.base,a.out)
