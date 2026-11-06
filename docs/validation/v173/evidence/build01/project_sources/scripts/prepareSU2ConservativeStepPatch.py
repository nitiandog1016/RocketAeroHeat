"""模块用途。在v164准确源上生成有限物性限步和只读事件上下文，所有终态恢复共享源。"""
from contextlib import contextmanager
from pathlib import Path
from buildSU2IdealAirRuntime import ROOT, SOURCE, BUILD, sha, save
from prepareSU2IdealAirPatch import replace
from su2ThermalRuntime import read

PARENT=ROOT/'outputs/v164_validation/build01'
FILE='SU2_CFD/src/variables/CNSVariable.cpp'
SOLVER='SU2_CFD/src/solvers/CNSSolver.cpp'
HEADER='SU2_CFD/include/variables/ConservativeStepTrace.hpp'
FLAVOR='rah-ideal-air-step-v1'


def patch_variable(text):
    text=replace(text,'#include "../../include/variables/AdmissibilityTrace.hpp"',
                 '#include "../../include/variables/ConservativeStepTrace.hpp"')
    start=text.index('bool CNSVariable::SetPrimVar(')
    tail=text.index('  /*--- Set enthalpy ---*/',start)
    return text[:start]+(ROOT/'third_party/ideal_air/ConservativeStep.inc').read_text('utf8')+text[tail:]


def patch_solver(text):
    text=replace(text,'#include "../../include/variables/CNSVariable.hpp"',
        '#include "../../include/variables/CNSVariable.hpp"\n#include "../../include/variables/ConservativeStepTrace.hpp"')
    marker='  /*--- Common preprocessing steps (implemented by CEulerSolver) ---*/'
    context='''  // 模块用途。只在事件发生时读真实几何/调用阶段/现存线性量，不改变更新。
  rah_step::ContextGuard step_context([&, this](unsigned long point) {
    rah_step::Context c;
    c.known=true; c.domain=point<nPointDomain; c.output=Output;
    c.global=geometry->nodes->GetGlobalIndex(point);
    c.inner=config->GetInnerIter(); c.time=config->GetTimeIter(); c.mesh=iMesh; c.rk=iRKStep;
    c.volume=geometry->nodes->GetVolume(point);
    c.alpha=nodes->GetUnderRelaxation(point);
    for(unsigned long j=0;j<nDim;++j) c.coordinates[j]=geometry->nodes->GetCoord(point,j);
    for(unsigned long j=0;j<nVar;++j) {
      const auto column=(j==nVar-1)?4:j;
      c.linear[column]=LinSysSol(point,j);
      c.residual[column]=LinSysRes(point,j);
    }
    return c;
  });

'''
    return replace(text,marker,context+marker)


@contextmanager
def overlay(folder):
    old_tree=read(PARENT/'source_tree.json')
    snapshots={SOURCE/p.relative_to(PARENT/'executed_sources'):p.read_bytes()
               for p in (PARENT/'executed_sources').rglob('*') if p.is_file()}
    snapshots[SOURCE/FILE]=patch_variable(snapshots[SOURCE/FILE].decode('utf8')).encode('utf8')
    snapshots[SOURCE/SOLVER]=patch_solver(snapshots[SOURCE/SOLVER].decode('utf8')).encode('utf8')
    snapshots[SOURCE/HEADER]=(ROOT/'third_party/ideal_air/ConservativeStepTrace.hpp').read_bytes()
    originals={p:p.read_bytes() if p.exists() else None for p in snapshots}
    try:
        for path,raw in snapshots.items(): path.write_bytes(raw)
        assert all(sha(Path(n))==s for n,s in old_tree.items() if n not in (str(SOURCE/FILE),str(SOURCE/SOLVER)))
        tree={n:sha(Path(n)) for n in old_tree}
        tree[str(SOURCE/HEADER)]=sha(SOURCE/HEADER)
        save(folder/'source_tree.json',tree)
        for name in ('docs/validation/v165/PLAN.md','scripts/prepareSU2ConservativeStepPatch.py',
                     'third_party/ideal_air/ConservativeStepTrace.hpp','third_party/ideal_air/ConservativeStep.inc'):
            target=folder/'project_sources'/name
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes((ROOT/name).read_bytes())
        yield snapshots
        assert all(sha(Path(n))==s for n,s in tree.items())
    finally:
        restored={}
        for path,raw in originals.items():
            target=folder/'executed_sources'/path.relative_to(SOURCE)
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(path.read_bytes())
            if raw is None:
                path.resolve().relative_to(SOURCE.resolve())
                path.unlink()
                restored[str(path)]=not path.exists()
            else:
                path.write_bytes(raw)
                restored[str(path)]=path.read_bytes()==raw
        save(folder/'source_restore.json',restored)
        assert all(restored.values())
