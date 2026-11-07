"""v198原生最小检查：物性查询小工具、polyMesh及独立同EOS参考。"""
from pathlib import Path
import argparse,json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,BASE,CASE,sha,save
from openFoamAirThermo import mixture_dictionary,evaluate,contract,R,HF,absolute_properties


def native():
    """模块一：隔离构建物性查询器，不编译引擎；各项均保留准确源和回执。"""
    from runOpenFoamProcess import run
    start=time.monotonic();out=BASE/'native_entry01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    build=out/'probe';shutil.copytree(ROOT/'native/openfoamThermoProbe',build)
    (build/'mixture.dict').write_text(mixture_dictionary())
    locks=[Path(__file__),ROOT/'scripts/openFoamAirThermo.py',ROOT/'docs/validation/v198/MINIMAL_ENTRY_PLAN.md']
    b=run(['wmake'],build,out/'build',90,locks+[build/'thermoProbe.C',build/'Make/files',build/'Make/options'],32)
    assert b['passed'],'物性查询器构建失败'
    q=run([build/'rahThermoProbe',build/'mixture.dict'],build,out/'query',20,
        locks+[build/'rahThermoProbe',build/'mixture.dict'],32)
    assert q['passed'],'原生物性查询失败'
    text=(out/'query/stdout.txt').read_text();assert 'out of temperature range' not in text
    rows=np.array([[float(x) for x in line.split()[1:]] for line in text.splitlines() if line.startswith('ROW ')])
    meta=next(line.split() for line in text.splitlines() if line.startswith('R '))
    assert len(rows)==9 and abs(float(meta[1])/R-1)<2e-15 and abs(float(meta[3])-HF)<2e-9
    errors={}
    fields=['cpJkgK','cvJkgK','enthalpyJkg','sensibleEnthalpyJkg','internalEnergyJkg','densityKgm3','soundSpeedMps']
    g=evaluate(5000,rows[:,0])
    for i,key in enumerate(fields):
        err=float(np.max(abs(rows[:,i+1]-g[key])/np.maximum(1,abs(g[key]))));errors[key]=err
        assert err<2e-10,key
    inv=float(np.max(abs(rows[:,-1]-rows[:,0])));assert inv<2e-7
    mesh=run(['checkMesh','-allTopology','-allGeometry'],CASE,out/'mesh',60,
        locks+[p for p in (CASE/'constant/polyMesh').iterdir()],64)
    meshlog=(out/'mesh/stdout.txt').read_text()
    assert mesh['passed'] and 'Mesh OK.' in meshlog,'原生网格审核失败'
    report=dict(passed=True,build=b,query=q,mesh=mesh,propertyMaxRelativeErrors=errors,
        inverseTemperatureMaxAbsK=inv,contract=contract(),elapsedSeconds=time.monotonic()-start,
        nativeRows=rows.tolist(),physicalQualification=False,CFDRuns=0,engineBuilds=0,helperBuilds=1)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','elapsedSeconds','propertyMaxRelativeErrors','inverseTemperatureMaxAbsK']}))


def independent():
    """模块二：导数/拒绝与守恒参照各有独立意义，不只复述配置值。"""
    from hybridEulerReference import control_reference
    start=time.monotonic();out=BASE/'independent01';out.mkdir(exist_ok=False)
    rejects=[]
    for p,t in [(0,300),(-1,300),(5000,199.99),(5000,2500.01),(5000,float('nan')),(float('inf'),300)]:
        try:evaluate(p,t)
        except ValueError:rejects.append([str(p),str(t)])
        else:raise AssertionError('超域未拒绝')
    ts=np.array([220.,300.,800.,1200.,2000.,2400.]);eps=.001
    cp=absolute_properties(ts)[0];der=(absolute_properties(ts+eps)[1]-absolute_properties(ts-eps)[1])/(2*eps)
    error=float(np.max(abs(der/cp-1)));assert error<2e-9
    ref=control_reference(eos=evaluate,gas_constant=R,model='rah-nasa7-air-of13-sensible-e-v1')
    assert max(v for k,v in ref['residuals'].items() if k.endswith('Relative'))<1e-10
    assert ref['residuals']['cpIntegralEnthalpyAbsolute']<1e-7
    # 拒绝不是把旧控制的网格资格转给新格式，实际checkMesh另验。
    from prepareOpenFoamControl import read_tetra,topology,OLD
    p,c,markers=read_tetra(OLD/'mesh.su2')
    missing={k:list(v) for k,v in markers.items()};missing['wall']=missing['wall'][:-1]
    duplicate={k:list(v) for k,v in markers.items()};duplicate['wall']=duplicate['wall']+[duplicate['wall'][0]]
    bad=c.copy();bad[0]=bad[1]
    mesh_rejects=[]
    for name,cc,mm in [('missing-boundary',c,missing),('duplicate-boundary',c,duplicate),('duplicate-cell',bad,markers)]:
        try:topology(p,cc,mm)
        except AssertionError:mesh_rejects.append(name)
        else:raise AssertionError(name+'未拒绝')
    save(out/'reference.json',ref)
    result=dict(passed=True,thermoRejects=rejects,meshRejects=mesh_rejects,
        cpDerivativeMaxRelativeError=error,reference=ref,elapsedSeconds=time.monotonic()-start,
        physicalQualification=False,CFDRuns=0)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','elapsedSeconds','cpDerivativeMaxRelativeError']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['native','independent']);args=parser.parse_args()
    (native if args.stage=='native' else independent)()
