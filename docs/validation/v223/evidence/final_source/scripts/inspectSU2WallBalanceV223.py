"""模块职责：一次读取准确热壁链及四份已有壁节点，形成完整接受决定的证据。"""
from pathlib import Path
import csv,json,shutil,sys,time
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,sha,read,save
BASE=ROOT/'outputs/v223_validation';DOC=ROOT/'docs/validation/v223'
OLD=ROOT/'tools/runtime/su2-thermal-build/source';NEW=ROOT/'outputs/v220_validation/source'

def main():
    out=BASE/'inspection01';out.mkdir(exist_ok=False);start=time.monotonic()
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    names=['SU2_CFD/src/integration/CIntegration.cpp','SU2_CFD/include/integration/CIntegration.hpp',
        'SU2_CFD/include/variables/CVariable.hpp','SU2_CFD/include/variables/CEulerVariable.hpp',
        'SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp','SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',
        'SU2_CFD/src/solvers/CNSSolver.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp',
        'SU2_CFD/src/variables/CNSVariable.cpp','Common/src/CConfig.cpp']
    sources=[]
    for name in names:
        p=NEW/name
        if p.is_file() and str(p) in deps:expected=deps[str(p)];assert sha(p)==expected
        else:
            expected=tree[str(OLD/name)]
            candidates=[ROOT/'outputs'/version/name for version in ('v222_validation/source01','v222_validation/source02',
                'v221_validation/source01','v219_validation/source01/source','v167_validation/build01/executed_sources')]+[OLD/name]
            p=next((q for q in candidates if q.is_file() and sha(q)==expected),None);assert p,name
        dest=out/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        sources.append(dict(relative=name,path=str(p),sha256=expected))
    g=np.load(str(ROOT/'outputs/v222_validation/geometry01/cache.npz'));ids=g['wallIDs'];tri=g['triangles'];area=g['faceArea']
    assert np.array_equal(ids,np.arange(len(ids))), 'This bounded wall-only reader requires existing contiguous original wall IDs'
    header=None;records=[];locks={}
    for step in (50,100,150,199):
        path=ROOT/'outputs/v222_validation/rocket01'/('restart_flow_%06d.csv'%step)
        expected=read(ROOT/'outputs/v222_validation/audit01/report.json')['outputSHA256'][path.name]
        assert sha(path)==expected;locks[str(path)]=expected
        with path.open(encoding='utf8') as stream:
            reader=csv.reader(stream,skipinitialspace=True);names=next(reader);ix={k:i for i,k in enumerate(names)}
            keys=['PointID','Density','Momentum_x','Momentum_y','Momentum_z','Energy','Pressure','Temperature',
                'Nu_Tilde','Residual_Density','Residual_Momentum_x','Residual_Momentum_y','Residual_Momentum_z','Residual_Energy','Residual_Nu_Tilde']
            a=np.array([[float(row[ix[k]]) for k in keys] for _,row in zip(range(len(ids)),reader)])
        assert np.array_equal(a[:,0],ids) and np.isfinite(a).all();v={k:a[:,j] for j,k in enumerate(keys)}
        speed=np.sqrt(sum(v['Momentum_'+d]**2 for d in ('x','y','z')))/v['Density']
        dt=v['Temperature']-300;face_dt=dt[tri].mean(axis=1)
        rms=np.sqrt(np.mean(np.column_stack([v['Residual_'+n] for n in ('Density','Momentum_x','Momentum_y','Momentum_z','Energy')])**2,axis=0))
        re=v['Residual_Energy']-300000*v['Residual_Density']
        row=dict(iteration=step,wallNodes=len(ids),wallSpeedMaxMps=float(max(speed)),wallNuTildeMaxAbs=float(max(abs(v['Nu_Tilde']))),
            wallTemperatureDeviationRangeK=[float(min(dt)),float(max(dt))],wallNodeTemperatureRmsDeviationK=float(np.sqrt(np.mean(dt**2))),
            originalFaceTemperatureAreaRmsDeviationK=float(np.sqrt(area@(face_dt**2)/sum(area))),
            faceAreaAbsMeanTemperatureErrorAbove1K_M2=float(sum(area[abs(face_dt)>1])),
            wallNativeRHS_RMS=rms.tolist(),wallPhysicalEnergyRHS_RMS=float(np.sqrt(np.mean(re**2))),
            noSlipOnSavedState=bool(max(speed)<1e-12),saWallZero=bool(max(abs(v['Nu_Tilde']))==0),physicalQualification=False)
        records.append(row)
        save(out/('wall_%06d.json'%step),row)
    assert all(sha(n)==h for n,h in locks.items())
    report=dict(passed=all(r['noSlipOnSavedState'] and r['saWallZero'] for r in records),sources=sources,records=records,
        sourceFrameSHA256=locks,sourceGeometrySHA256=sha(ROOT/'outputs/v222_validation/geometry01/cache.npz'),
        strongTemperatureDirichlet=False,weakWallTemperatureUsesGiven300K=True,fullBoundaryBalanceObserved=False,
        savedWallVelocityFailure=False,uniqueCauseEstablished=False,CFDRuns=0,nativeQueries=0,physicalQualification=False,
        seconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(dict(passed=report['passed'],wallSpeedMax=max(r['wallSpeedMaxMps'] for r in records),seconds=report['seconds'])))

def localize():
    """仅读尚未保存的逐点质量/能量RHS分布，判断后段增长是否主要来自热壁。"""
    out=BASE/'localize01';out.mkdir(exist_ok=False)
    audit=read(ROOT/'outputs/v222_validation/audit01/report.json');g=np.load(str(ROOT/'outputs/v222_validation/geometry01/cache.npz'))
    nwall=len(g['wallIDs']);npoint=len(g['points']);records=[]
    for step in (50,100,150,199):
        p=ROOT/'outputs/v222_validation/rocket01'/('restart_flow_%06d.csv'%step)
        assert sha(p)==audit['outputSHA256'][p.name]
        with p.open(encoding='utf8') as stream:header=next(csv.reader(stream,skipinitialspace=True))
        keys=['PointID','x','y','z','Residual_Density','Residual_Energy']
        a=np.loadtxt(str(p),delimiter=',',skiprows=1,usecols=[header.index(k) for k in keys])
        assert np.array_equal(a[:,0],np.arange(npoint)) and np.isfinite(a).all()
        mass=a[:,4];energy=a[:,5]-300000*mass
        fraction=lambda x:float(np.sum(x[:nwall]**2)/np.sum(x**2))
        rms=lambda x:float(np.sqrt(np.mean(x*x)))
        top=np.argsort(abs(energy))[-20:][::-1]
        row=dict(iteration=step,wallMassSquaredResidualFraction=fraction(mass),wallEnergySquaredResidualFraction=fraction(energy),
            totalMassRMS=rms(mass),wallMassRMS=rms(mass[:nwall]),nonwallMassRMS=rms(mass[nwall:]),
            totalPhysicalEnergyRMS=rms(energy),wallPhysicalEnergyRMS=rms(energy[:nwall]),nonwallPhysicalEnergyRMS=rms(energy[nwall:]),
            largestEnergyNodes=[dict(id=int(i),xyz=a[i,1:4].tolist(),physicalEnergyRHS=float(energy[i]),massRHS=float(mass[i]),wall=bool(i<nwall)) for i in top],
            stage='pre-update implicit RHS, not same-stage saved primitive or complete momentum balance',uniqueCauseEstablished=False)
        records.append(row);assert sha(p)==audit['outputSHA256'][p.name]
    save(out/'report.json',dict(passed=True,records=records,sourceAuditSHA256=sha(ROOT/'outputs/v222_validation/audit01/report.json'),
        physicalQualification=False,CFDRuns=0,nativeQueries=0))
    print(json.dumps(dict(passed=True,lastWallEnergyFraction=records[-1]['wallEnergySquaredResidualFraction'])))

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    elif len(sys.argv)>1 and sys.argv[1]=='task_localize':localize()
    else:
        # 复用已测进程监督，仅改变本轮证据目录；整个本轮输出作终态独立容量检查。
        BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
        real_size=guard.size_bytes;guard.size_bytes=lambda unused:real_size(BASE)
        locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            ROOT/'outputs/v222_validation/audit01/report.json',ROOT/'outputs/v222_validation/geometry01/cache.npz']
        stage='localize' if len(sys.argv)>1 and sys.argv[1]=='localize' else 'inspect'
        raise SystemExit(0 if guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_localize' if stage=='localize' else 'task'],120,128,locks,BASE,'offline') else 1)
