"""依据成熟限制器源码及初态面耗散证据，准备唯一同场空间离散对照。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from phaseAwareViscousAir import PhaseAwareViscousAir

def main():
    start=time.monotonic();base=ROOT/'outputs/v207_validation';out=base/'limited_prepare02';out.mkdir(exist_ok=False)
    old=base/'thermal02';case=base/'thermal03'
    if case.exists() and any(case.iterdir()):raise ValueError('不能覆盖已有候选')
    case.mkdir(exist_ok=True)
    first=json.loads((base/'thermal_audit01/report.json').read_text())
    if first['passed'] or not first['nativeSignAndGradientInterfacePassed']:raise ValueError('不符已识别的空间离散缺项')
    engine=Path('/opt/openfoam13');names=['src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/vanLeer/vanLeer.H',
        'src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/NVDTVD.H',
        'tutorials/shockFluid/wedge15Ma5/system/fvSchemes','tutorials/shockFluid/wedge15Ma5/system/controlDict']
    sources=[]
    for name in names:
        p=engine/name;q=out/'source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        sources.append(dict(path=str(p),sha256=sha(p)))
    limiter=(engine/names[0]).read_text();nvd=(engine/names[1]).read_text()
    if '(r + mag(r))/(1 + mag(r))' not in limiter or '2*(gradcf/gradf) - 1' not in nvd:raise ValueError('限制器公式不同')
    for d in ['0','constant','system']:shutil.copytree(old/d,case/d)
    shutil.copyfile(old/'thermo_contract.json',case/'thermo_contract.json')
    path=case/'system/fvSchemes';s=path.read_text()
    for a,b in [('reconstruct(rho) upwind','reconstruct(rho) vanLeer'),('reconstruct(U) upwind','reconstruct(U) vanLeerV'),('reconstruct(T) upwind','reconstruct(T) vanLeer')]:
        if s.count(a)!=1:raise ValueError('原重构配置改变')
        s=s.replace(a,b)
    path.write_text(s)
    unchanged={str(p.relative_to(old)):sha(p)==sha(case/p.relative_to(old)) for d in ['0','constant','system'] for p in (old/d).rglob('*') if p.is_file() and p.relative_to(old).as_posix()!='system/fvSchemes'}
    if not all(unchanged.values()):raise ValueError('指定空间重构以外发生变化')
    gas=PhaseAwareViscousAir();oracle=json.loads((base/'period_repair01/oracle.json').read_text());T=np.array(oracle['temperatureK'])
    rho=gas.evaluate(100000,T)['densityKgm3'];dy=.001/32
    ends=gas.evaluate(100000,np.array([300.,600.]))['densityKgm3']
    faces=np.r_[ends[0],.5*(rho[:-1]+rho[1:]),ends[1]];grad=np.diff(faces)/dy;delta=np.diff(rho)
    lim=lambda r:(r+abs(r))/(1+abs(r))
    left=rho[:-1]+.5*lim(2*grad[:-1]*dy/delta-1)*delta
    right=rho[1:]-.5*lim(2*grad[1:]*dy/delta-1)*delta
    ratio=float(abs(left-right).max()/abs(delta).max())
    if ratio>.01:raise ValueError('事前耗散降低预测未过')
    shutil.copyfile(__file__,out/'executed_entry.py')
    report=dict(passed=True,sources=sources,maxDensityFaceJumpRatio=ratio,unchangedInputs=unchanged,
        sourceModel='official Euler/vanLeer/vanLeerV/Kurganov with stock implicit viscous terms',
        inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in case.rglob('*') if p.is_file()},
        originalHeatControlPassed=False,CFDRuns=0,builds=0,meshReruns=0,elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    save(out/'report.json',report);print(json.dumps(dict(passed=True,densityJumpRatio=ratio,seconds=report['elapsedSeconds'])))

if __name__=='__main__':main()
