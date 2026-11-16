"""模块职责：在已有四管实际经过区域分离源场和插值误差，并绑定准确成熟原生源码。"""
from pathlib import Path
import json,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh
from auditForebodyStreamtubesV241 import independent_states
BASE=ROOT/'outputs/v243_validation';DOC=ROOT/'docs/validation/v243';OLD=ROOT/'outputs/v242_validation'
CLOSED=True


def inspect():
    out=BASE/'inspect02';out.mkdir(exist_ok=False)
    manifest=read(OLD/'tubes01/manifest.json');audit=read(OLD/'flow_audit01/report.json');assert read(OLD/'tube_audit01/report.json')['passed']
    x,t,markers=read_tetra_mesh(OLD/'run01/mesh.su2');raw=np.loadtxt(str(OLD/'run01/state_000400.csv'),delimiter=',',skiprows=1);raw=raw[np.argsort(raw[:,0])]
    assert np.array_equal(raw[:,1:4],x);u=raw[:,4:9];records=[];cellsets=[];stopfaces=[];max_interpolation=0.
    # 模块一：只访问实际射线经过单元的供体节点；两种H取法均记录，不替换生产状态。
    for info in manifest['tubes']:
        z=np.load(str(OLD/'tubes01'/info['file']));ids=np.unique(np.concatenate([z['ray%d_cell'%j] for j in range(9)]));cellsets.append(ids)
        donor=np.unique(t[ids]);states=independent_states(u[donor]);hn=np.full(len(x),np.nan);hn[donor]=states['H'];sn=np.full(len(x),np.nan);sn[donor]=states['S']
        rows=[]
        for j in range(9):
            cells=z['ray%d_cell'%j];weights=z['ray%d_weights'%j];sample=independent_states(z['ray%d_U'%j]);n=len(cells)
            linearH=np.sum(weights*hn[t[cells]],axis=1);linearS=np.sum(weights*sn[t[cells]],axis=1)
            nonlinearH=sample['H']-linearH;nonlinearS=sample['S']-linearS
            max_interpolation=max(max_interpolation,float(max(abs(nonlinearH))/abs(audit['freestreamHJkg'])))
            rows.append(np.column_stack([np.full(n,j),np.arange(n),cells,z['ray%d_s'%j],sample['H'],linearH,nonlinearH,sample['S'],linearS,nonlinearS]))
            stop=info['rayStops'][j]
            if stop['stopReason']=='TUBE_BOUNDARY_CROSSING' and tuple(sorted(stop['stopBoundaryNodes'])) in set(tuple(sorted(f)) for f in markers['wall']):stopfaces.append(stop['stopBoundaryNodes'])
        data=np.vstack(rows);np.savetxt(str(out/('tube%02d_decomposition.csv'%info['id'])),data,delimiter=',',fmt='%.17g',comments='',header='ray,station,cell,sM,Hsample,Hlinear,Hmixing,Ssample,Slinear,Smixing')
        center=rows[0][:info['commonStations']];den=abs(center[0,4]);dtotal=center[:,4]-center[0,4];dsource=center[:,5]-center[0,5];dmix=center[:,6]-center[0,6]
        assert np.max(abs(dtotal-dsource-dmix))<2e-9
        records.append(dict(id=info['id'],cells=len(ids),donorNodes=len(donor),commonStations=len(center),
          sourceNodeHRelativeRange=[float(hn[donor].min()/audit['freestreamHJkg']-1),float(hn[donor].max()/audit['freestreamHJkg']-1)],
          centerMaxTotalHDrift=float(max(abs(dtotal))/den),centerMaxLinearNodeHDrift=float(max(abs(dsource))/den),
          centerMaxMixingHDrift=float(max(abs(dmix))/den),decompositionClosureMaxJkg=float(max(abs(dtotal-dsource-dmix)))))
    # 模块二：P1质量向量的真实常散度和跨壁通量；不是SU2原生数值面通量。
    ids=np.unique(np.concatenate(cellsets));p=x[t[ids]];inverse=np.linalg.inv((p[:,1:]-p[:,0,None]).transpose(0,2,1))
    gradients=np.concatenate([-inverse.sum(axis=1)[:,None,:],inverse],axis=1)
    momentum=u[t[ids],1:4];divergence=np.einsum('nkd,nkd->n',momentum,gradients)
    volume=abs(np.linalg.det((p[:,1:]-p[:,0,None])))/6;speed=np.linalg.norm(momentum.mean(axis=1),axis=1)
    scaled=abs(divergence)*np.cbrt(6*volume)/np.maximum(speed,1e-30)
    faces=np.unique(np.sort(np.asarray(stopfaces,np.int64),axis=1),axis=0);q=x[faces]
    normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])/2;mass=u[faces,1:4].mean(axis=1)
    through=np.einsum('ij,ij->i',normal,mass);fraction=abs(through)/(np.linalg.norm(normal,axis=1)*np.linalg.norm(mass,axis=1))
    np.savez_compressed(str(out/'transport.npz'),cellIDs=ids,divMass=divergence,scaledDivMass=scaled,wallFaces=faces,p1NormalMassFlow=through,p1NormalFraction=fraction)
    # 模块三：读取与真实驱动依赖逐SHA相同的准确源，不拿当前其他分支替代。
    dependencies=read(ROOT/'outputs/v229_validation/process_recover_build/report.json')['inputSHA256'];source=ROOT/'outputs/v220_validation/source/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl'
    assert sha(source)==dependencies[str(source)]
    text=source.read_text(encoding='utf8').splitlines();lo=next(i for i,line in enumerate(text) if '::BC_Sym_Plane(' in line);hi=next(i for i in range(lo+1,len(text)) if '::BC_Periodic(' in text[i])
    (out/'native_slip_excerpt.txt').write_text('\n'.join('%d: %s'%(i+1,text[i]) for i in range(lo,hi))+'\n',encoding='utf8')
    report=dict(passed=True,sourceUpdate=400,actualTubes=records,visitedUniqueCells=len(ids),visitedCellMaxScaledDivergence=float(max(scaled)),
      visitedCellScaledDivergencePercentiles=np.percentile(scaled,[0,50,95,100]).tolist(),actualStoppedWallFaces=len(faces),
      p1WallNormalMassFractionRange=[float(min(fraction)),float(max(fraction))],maxLocalHInterpolationOffsetRelativeHInf=max_interpolation,
      nativeSlipSource=str(source),nativeSlipSHA256=sha(source),nativeSlipLines=[lo+1,hi],
      P1DivergenceIsNotNativeResidual=True,wallFluxIsNotNativeRoeFlux=True,sourceValuesChanged=False,CFDUpdates=0,
      trajectoryReintegrations=0,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'report.json',report);print(json.dumps(report))


def main():
    if CLOSED:raise SystemExit('v243 fixed inspection closed; use NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
    inputs=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/auditForebodyStreamtubesV241.py',
      ROOT/'data/thermo/phase_aware_dry_air_v1.json',OLD/'flow_audit01/report.json',OLD/'tube_audit01/report.json',OLD/'flow_prepare01/dependencies.json',
      ROOT/'outputs/v220_validation/source/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',ROOT/'outputs/v229_validation/process_recover_build/report.json']+[OLD/'run01'/n for n in ('mesh.su2','state_000400.csv')]+list((OLD/'tubes01').glob('*'))
    return guard.supervise('inspect2',[sys.executable,str(Path(__file__).resolve()),'task'],60,256,inputs,BASE,'offline')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v243 fixed inspection closed; use NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':inspect()
    else:raise SystemExit(0 if main() else 1)
