"""模块用途。只在准确admissible-v1的SST共同松弛入口增加热内能约束。"""
from contextlib import contextmanager
from pathlib import Path
from buildSU2IdealAirRuntime import ROOT,SOURCE,BUILD
from prepareSU2IdealAirPatch import replace
from su2ThermalRuntime import read,sha
from validateIdealAir import save

PARENT=ROOT/'outputs/v164_validation/build01'
FILE='SU2_CFD/src/solvers/CTurbSSTSolver.cpp'
HEADER='SU2_CFD/include/fluid/SSTEnergyUpdate.hpp'
FLAVOR='rah-sst-energy-v1'


def patch(text):
    text=replace(text,'#include "../../include/solvers/CTurbSSTSolver.hpp"',
        '#include "../../include/solvers/CTurbSSTSolver.hpp"\n#include "../../include/fluid/SSTEnergyUpdate.hpp"')
    needle='  ComputeUnderRelaxationFactorHelper(solver_container, allowableRatio);\n}'
    extra='''  ComputeUnderRelaxationFactorHelper(solver_container, allowableRatio);
  // 模块用途。原SST共同松弛之后，按真实旧解/当前密度/线性增量约束热内能。
  if (config->GetKind_FluidModel() != IDEAL_AIR_VARIABLE_CP) return;
  const auto flowNodes = solver_container[FLOW_SOL]->GetNodes();
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint=0; iPoint<nPointDomain; ++iPoint) {
    rah_sst_energy::Input input;
    input.dimension=nDim; input.flow[0]=flowNodes->GetSolution(iPoint,0);
    for (unsigned j=0;j<nDim;++j) input.flow[j+1]=flowNodes->GetSolution(iPoint,j+1);
    input.flow[4]=flowNodes->GetSolution(iPoint,nDim+1);
    input.native_alpha=nodes->GetUnderRelaxation(iPoint);
    for (unsigned j=0;j<2;++j) {
      input.old[j]=nodes->GetSolution_Old(iPoint,j); input.increment[j]=LinSysSol(iPoint,j);
      input.lower[j]=lowerlimit[j]; input.upper[j]=upperlimit[j];
    }
    rah_sst_energy::Result result;
    try {
      result=rah_sst_energy::constrain(input);
    } catch (const std::domain_error&) {
      rah_sst_energy::record(iPoint,input,result,"invalid-base-or-input");
      throw;
    }
    if(result.limited) rah_sst_energy::record(iPoint,input,result,"limited");
    nodes->SetUnderRelaxation(iPoint,result.alpha);
  }
  END_SU2_OMP_FOR
}'''
    return replace(text,needle,extra)


@contextmanager
def overlay(folder):
    tree=read(PARENT/'source_tree.json')
    snapshots={SOURCE/p.relative_to(PARENT/'executed_sources'):p.read_bytes()
               for p in (PARENT/'executed_sources').rglob('*') if p.is_file()}
    assert sha(SOURCE/FILE)==tree[str(SOURCE/FILE)]
    snapshots[SOURCE/FILE]=patch((SOURCE/FILE).read_text('utf8')).encode('utf8')
    snapshots[SOURCE/HEADER]=(ROOT/'third_party/ideal_air/SSTEnergyUpdate.hpp').read_bytes()
    originals={p:p.read_bytes() if p.exists() else None for p in snapshots}
    try:
        for p,raw in snapshots.items():p.write_bytes(raw)
        assert all(sha(Path(n))==digest for n,digest in tree.items() if n!=str(SOURCE/FILE))
        current={n:sha(Path(n)) for n in tree};current[str(SOURCE/HEADER)]=sha(SOURCE/HEADER)
        save(folder/'source_tree.json',current)
        yield snapshots
        assert all(sha(Path(n))==digest for n,digest in current.items())
    finally:
        restored={}
        for p,raw in originals.items():
            target=folder/'executed_sources'/p.relative_to(SOURCE);target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(p.read_bytes())
            if raw is None:p.resolve().relative_to(SOURCE.resolve());p.unlink();restored[str(p)]=not p.exists()
            else:p.write_bytes(raw);restored[str(p)]=p.read_bytes()==raw
        save(folder/'source_restore.json',restored);assert all(restored.values())
