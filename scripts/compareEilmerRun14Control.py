"""将具名Run14原生载荷映射至原测站，输出有符号全周功率与未收敛对照。"""
from pathlib import Path
import argparse
import csv
import json
import math
import re
import time
import numpy as np
from readEilmerRun14Control import grid,field,need
from runEilmerProcess import sha


def tangential_wall_shear(flux, outsign, tangent, density, viscosity, distance):
    """原生列是完整黏性动量通量；投影切向后计算壁面摩擦速度，排除法向应力。"""
    need(all(np.isfinite(z).all() for z in (flux,outsign,tangent,density,viscosity,distance)), '剪切输入非有限')
    need(np.all(density>0) and np.all(viscosity>0) and np.all(distance>0), '壁物性/法向距离非正')
    signed=outsign*np.sum(flux*tangent,axis=1)
    return signed,distance*np.sqrt(abs(signed)/density)/(viscosity/density)


def interpolate_station(x,value,block,at):
    """同几何段内插值；拐角/端点外或数据非有限一律无值，不用0填充。"""
    need(np.isfinite(at),'测站坐标非有限')
    target=0 if at<.1017 else 1
    take=block==target;xx=x[take];vv=value[take]
    need(len(xx)>1 and np.all(np.diff(xx)>0),'源曲线重复/乱序')
    if at<xx[0] or at>xx[-1]:return None,'outside-computed-face-centre-support'
    j=int(np.searchsorted(xx,at))
    if j<len(xx) and at==xx[j]:
        return (float(vv[j]),'matched-centre') if np.isfinite(vv[j]) else (None,'unknown-value')
    if not np.isfinite(vv[j-1:j+1]).all():return None,'unknown-bracketing-value'
    return float(vv[j-1]+(at-xx[j-1])/(xx[j]-xx[j-1])*(vv[j]-vv[j-1])),'linear-within-one-geometric-region'


def read_wall(case,index):
    curves=[]
    for b in (0,1):
        matches=list((Path(case)/('lmrsim/loads/%04d'%index)).glob('blk-%04d-*-wall.dat'%b))
        need(len(matches)==1,'缺少完整壁面载荷')
        a=np.genfromtxt(matches[0],names=True,deletechars='');need(len(a)==160,'壁面数量错误')
        need(all(np.isfinite(a[k]).all() for k in a.dtype.names),'非有限壁面数据')
        g=grid(case,b);p=g[0,:,:2];edge=np.diff(p,axis=0);ds=np.linalg.norm(edge,axis=1)
        center=(p[:-1]+p[1:])/2;normal=np.column_stack([-edge[:,1],edge[:,0]])/ds[:,None]
        pos=np.column_stack([a[k] for k in ['pos.x','pos.y']]);nn=np.column_stack([a[k] for k in ['n.x','n.y']])
        need(np.max(abs(center-pos))<3e-14,'原生面位置不符')
        # 轴对称内部面积为每弧度r*ds，必须乘2π才能得到完整周向面积/功率。
        need(np.max(abs(a['area']/(center[:,1]*ds)-1))<2e-11,'每弧度面积不符')
        need(np.max(abs(nn*a['outsign'][:,None]+normal))<2e-12,'原生法向/外向符号不符')
        need(np.max(abs(a['T']-295.2))<1e-10,'热壁温度错误')
        need(np.max(abs(a['q_total']-a['q_cond']-a['q_diff']))<2e-12*max(1.,float(np.max(abs(a['q_total'])))),'热分项不闭合')
        cells=(g[:-1,:-1,:2]+g[:-1,1:,:2]+g[1:,:-1,:2]+g[1:,1:,:2])/4
        distance=np.sum((cells[0]-center)*normal,axis=1)
        tau=np.column_stack([a[k] for k in ['tau_wall.x','tau_wall.y']])
        tangential,yplus=tangential_wall_shear(tau,a['outsign'],edge/ds[:,None],a['rho'],a['mu'],distance)
        curves.append(np.column_stack([center[:,0],center[:,1],np.full(160,b),2*math.pi*a['area'],a['p'],
            -a['outsign']*a['q_total'],a['T'],distance,yplus,tangential]))
    return np.concatenate(curves)


def compare(case,run_receipt,reference,out):
    start=time.monotonic();case=Path(case).resolve();out=Path(out);need(not out.exists(),'拒绝覆盖');out.mkdir(parents=True)
    locked=[Path(__file__),Path(__file__).with_name('readEilmerRun14Control.py'),Path(run_receipt),
            Path(reference)/'heat_flux_SI.csv',Path(reference)/'pressure_SI.csv']
    input_hashes={str(p.resolve()):sha(p) for p in locked}
    receipt=json.loads(Path(run_receipt).read_text());need(receipt['stage']=='run' and receipt['case']==str(case),'运行来源不符')
    need(receipt['caseManifestSHA256']==sha(case/'case_manifest.json'),'工况身份改变')
    for p,digest in receipt['outputSHA256'].items():need(sha(p)==digest,'保存结果变化')
    native_log=(Path(run_receipt).parent/'process00/stdout.txt').read_text(errors='replace')
    native_stop=re.findall(r'^STOP-REASON:\s*(.*)$',native_log,re.M)
    need(receipt['commands'] and receipt['commands'][0]['inputsStable'],'原任务输入变化')
    meta=np.atleast_2d(np.loadtxt(case/'lmrsim/loads/loads-metadata',skiprows=1));need(meta.shape[1]==3 and len(meta)>0,'无载荷帧')
    index=int(meta[-1,0]);step=int(meta[-1,2]);a=read_wall(case,index)
    headers=['x_m','r_m','block','full_ring_area_m2','pressure_Pa','heat_into_wall_W_m2','Tw_given_K','wall_normal_distance_m','independent_y_plus','signed_tangent_stress_Pa']
    with (out/'wall.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(a)
    measurements=[];metrics={}
    for quantity,filename,column in [('q','heat_flux_SI.csv',5),('p','pressure_SI.csv',4)]:
        errors=[];observed=[];predicted=[]
        for r in csv.DictReader((Path(reference)/filename).open()):
            if r['record_status']!='source-active':continue
            x=float(r['x_from_leading_edge_m']);value,reason=interpolate_station(a[:,0],a[:,column],a[:,2],x)
            truth=float(r['value_SI']) if r['value_SI'] else None
            if truth is None:value=None;reason='missing-measurement'
            error=None if value is None else value-truth
            measurements.append(dict(quantity=quantity,sensor=r['sensor'],x_m=x,measurement=truth,prediction=value,
                signedError=error,relativeError=None if error is None or truth==0 else error/truth,
                status=reason,experimentalUncertainty=None,physicsQualified=False,nativeStep=step))
            if error is not None:errors.append(error);observed.append(truth);predicted.append(value)
        obs=np.array(observed);err=np.array(errors);pred=np.array(predicted)
        metrics[quantity]=dict(matchedStations=len(obs),rmsAbsolute=float(np.sqrt(np.mean(err**2))),
            rmsDividedByMeasuredPeak=float(np.sqrt(np.mean(err**2))/max(abs(obs))),
            meanAbsoluteRelativeError=float(np.mean(abs(err/obs))),measuredPeakAtStations=float(obs.max()),
            predictedPeakAtStations=float(pred.max()),calculatedContinuousFacePeak=float(a[:,column].max()),
            calculatedPeakX=float(a[np.argmax(a[:,column]),0]))
    with (out/'stations.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(measurements[0]));w.writeheader();w.writerows(measurements)
    h=np.genfromtxt(case/'lmrsim/diagnostics/nk-diagnostics',names=True,deletechars='');h=np.atleast_1d(h)
    need(all(np.isfinite(h[k]).all() for k in h.dtype.names),'历史非有限')
    full=h[h['phase']==1];improvements={}
    if len(full)>=200:
        for key in ['mass-rel','total-energy-rel']:
            improvements[key]=float(np.mean(np.log10(full[key][:100]))-np.mean(np.log10(full[key][-100:])))
    stable={}
    if len(meta)>=2:
        old=read_wall(case,int(meta[-2,0]))
        stable['savedSteps']=[int(meta[-2,2]),step]
        need(np.array_equal(a[:,:4],old[:,:4]),'保存帧几何/面积改变')
        for q,j in [('p',4),('q',5)]:stable[q+'AreaL1RelativeChange']=float(np.dot(abs(a[:,j]-old[:,j]),a[:,3])/max(1e-300,np.dot(abs(a[:,j]),a[:,3])))
        stable['qPeakRelativeChange']=float(abs(a[:,5].max()-old[:,5].max())/max(1e-300,abs(a[:,5].max())))
    # 按实际保存槽读取，不假定槽号等于迭代数。EOS与能量参考只是代数检查。
    snapshots=[]
    restart=np.atleast_2d(np.loadtxt(case/'lmrsim/snapshots/restart-metadata'))
    for p in sorted((case/'lmrsim/snapshots').iterdir()):
        if p.is_dir() and p.name.isdigit():
            fields=[field(case,b,int(p.name)) for b in (0,1)]
            snapshots.append(dict(slot=int(p.name),Tmin=min(float(z['T'].min()) for z in fields),Tmax=max(float(z['T'].max()) for z in fields)))
    crossings=[]
    for j in range(len(a)-1):
        if a[j,2]==a[j+1,2] and a[j,9]*a[j+1,9]<0:
            crossings.append(dict(xBracketM=a[j:j+2,0].tolist(),
                direction='positive-to-negative' if a[j,9]>0 else 'negative-to-positive',
                status='unconverged-shear-sign-change-not-validated-separation-or-reattachment'))
    result=dict(nativeProcessPassed=receipt['passed'],lastSavedLoadsStep=step,historyLastStep=int(h['step'][-1]),
        nativeStopReason=native_stop[-1] if native_stop else None,
        supervisorStopReason=receipt['commands'][0].get('stopReason'),
        reconstructionExtremaClipping=json.loads((case/'lmrsim/config').read_text())['extrema_clipping'],
        comparisonInterfacePassed=True,physicsQualified=False,newQualifiedAreaM2=0,
        fullRingWallAreaM2=float(a[:,3].sum()),researchHeatPowerIntoWallW=float(np.dot(a[:,3],a[:,5])),
        positiveHeatFaces=int(np.sum(a[:,5]>0)),negativeHeatFaces=int(np.sum(a[:,5]<0)),
        independentYPlus99=float(np.quantile(a[:,8],.99)),negativeTangentialStressFaces=int(np.sum(a[:,9]<0)),
        lastGlobalRelativeResidual=float(h['global-residual-rel'][-1]),
        lastMassRelativeResidual=float(h['mass-rel'][-1]),lastEnergyRelativeResidual=float(h['total-energy-rel'][-1]),
        meetsOriginalResidualTarget=bool(max(h[k][-1] for k in ['global-residual-rel','mass-rel','total-energy-rel'])<=1e-8),
        sameOrderLogImprovements=improvements,developmentBenefit=bool(len(improvements)==2 and min(improvements.values())>=.5 and len(meta)>=2),
        savedWallStability=stable,metrics=metrics,snapshots=snapshots,
        restartMetadataSteps=restart[:,0].astype(int).tolist(),
        yPlusDefinition='True face-normal centre distance and tangential viscous traction; wall rho/mu; not a mesh-convergence qualification',
        shearSignChangeBrackets=crossings,exactZeroShearFaces=int(np.sum(a[:,9]==0)),
        inputSHA256=input_hashes,inputsStable=all(sha(p)==v for p,v in input_hashes.items()),
        unmatched=[r['sensor'] for r in measurements if r['prediction'] is None],
        runReceiptSHA256=sha(run_receipt),elapsedSeconds=time.monotonic()-start,
        limitation='Unconverged research unless separate physics gates pass; pointwise sensors with unknown footprint/uncertainty; no interpolation across corner.')
    need(result['inputsStable'],'后处理输入变化')
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--run-receipt',required=True)
    p.add_argument('--reference',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    print(json.dumps(compare(a.case,a.run_receipt,a.reference,a.out)))
