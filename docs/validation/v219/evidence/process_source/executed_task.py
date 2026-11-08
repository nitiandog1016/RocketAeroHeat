"""模块职责：绑定v167已有SU2实数构建的定向源码，供整套NS/RANS方案一次审查。"""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v219_validation'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
NAMES=['Common/src/CConfig.cpp','Common/include/option_structure.hpp',
 'SU2_CFD/include/fluid/CSU2IdealAir.hpp','SU2_CFD/include/fluid/IdealAir.hpp',
 'SU2_CFD/include/fluid/CFluidModel.hpp','SU2_CFD/src/fluid/CFluidModel.cpp',
 'SU2_CFD/src/solvers/CSolverFactory.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp',
 'SU2_CFD/src/solvers/CNSSolver.cpp','SU2_CFD/src/variables/CNSVariable.cpp',
 'SU2_CFD/src/variables/CEulerVariable.cpp','SU2_CFD/src/solvers/CTurbSASolver.cpp',
 'SU2_CFD/src/solvers/CTurbSolver.cpp','SU2_CFD/include/solvers/CScalarSolver.inl',
 'SU2_CFD/include/solvers/CTurbSASolver.hpp','SU2_CFD/src/variables/CTurbSAVariable.cpp',
 'SU2_CFD/include/numerics/turbulent/turb_sources.hpp','SU2_CFD/include/numerics/turbulent/turb_sa_edge_flux.hpp',
 'SU2_CFD/src/integration/CSingleGridIntegration.cpp','SU2_CFD/src/iteration/CFluidIteration.cpp',
 'SU2_CFD/src/output/CFlowOutput.cpp','COPYING','AUTHORS.md']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def main():
    start=time.monotonic();out=BASE/'source01';out.mkdir(exist_ok=False)
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    runtime=ROOT/'tools/runtime/su2/8.6.0-rah-sst-energy-v1/runtime.json';r=read(runtime)
    assert sha(ROOT/'outputs/v167_validation/build01/source_tree.json')==r['sourceTreeSHA256']
    assert r['sourceCommit']=='7048fd35f56d4bf52e586429d6746a36a10255f8'
    rows=[]
    for name in NAMES:
        expected=tree[str(SOURCE/name)]
        candidates=[ROOT/('outputs/'+v+'/executed_sources')/name for v in (
            'v167_validation/build01','v164_validation/build01','v163_validation/build03')]+[SOURCE/name]
        matches=[p for p in candidates if p.is_file() and sha(p)==expected]
        assert matches,'No exact source for '+name
        p=matches[0];q=out/'source'/name;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(p.read_bytes())
        rows.append(dict(relative=name,origin=p.relative_to(ROOT).as_posix(),sha256=expected,bytes=p.stat().st_size,
                         path=q.relative_to(ROOT).as_posix()))
    executable=Path(r['executable']);expected=next(x['sha256'] for x in r['files'] if x['path']=='bin/SU2_CFD.exe')
    assert sha(executable)==expected
    report=dict(passed=True,sourceCommit=r['sourceCommit'],runtimePath=runtime.relative_to(ROOT).as_posix(),
        runtimeSHA256=sha(runtime),executableSHA256=expected,sourceTreeSHA256=r['sourceTreeSHA256'],sources=rows,
        sourceCount=len(rows),sourceBytes=sum(x['bytes'] for x in rows),nativeQueries=0,builds=0,CFDRuns=0,
        newCombinationAccepted=False,seconds=time.monotonic()-start)
    (out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'))
    print(json.dumps({k:report[k] for k in ('passed','sourceCount','sourceBytes','seconds','newCombinationAccepted')}))

if __name__=='__main__':main()
