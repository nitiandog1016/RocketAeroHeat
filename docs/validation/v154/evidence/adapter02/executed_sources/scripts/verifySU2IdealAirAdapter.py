"""编译真实SU2基类适配并核对准确原函数，保留执行源和有界进程回执。"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import zipfile

import numpy as np
from validateIdealAir import run
from verifySU2AUSMUp2Flux import body

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'tools/runtime/su2-thermal-build'
SOURCE=BASE/'source'
BUILD=BASE/'build-v117'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf8')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--attempt',required=True)
    args=parser.parse_args()
    assert args.attempt.isdigit()
    folder=ROOT/'outputs/v154_validation'/('adapter'+args.attempt)
    folder.mkdir(parents=True,exist_ok=False)
    archive=BASE/'su2.zip'
    assert sha(archive)=='c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    with zipfile.ZipFile(str(archive)) as bundle:
        name=next(n for n in bundle.namelist() if n.endswith('/SU2_CFD/src/numerics/CNumerics.cpp'))
        source=bundle.read(name).decode('utf8')
        diffusion_name=next(n for n in bundle.namelist() if n.endswith('/SU2_CFD/src/numerics/flow/flow_diffusion.cpp'))
        diffusion=bundle.read(diffusion_name).decode('utf8')
    signature='void CNumerics::GetInviscidProjJac(const su2double *val_velocity, const su2double *val_enthalpy,'
    fragment=body(source,signature)
    header='''#pragma once
struct NativeGeneralJacobian {
  unsigned short nDim=3;
  void GetInviscidProjJac(const double* val_velocity,const double* val_enthalpy,
      const double* val_chi,const double* val_kappa,const double* val_normal,
      double val_scale,double** val_Proj_Jac_Tensor) const
'''+fragment+'\n};\n'
    heat=body(diffusion,'void CGeneralAvgGrad_Flow::SetHeatFluxJacobian(')
    header+='''struct NativeGeneralHeatJacobian {
  unsigned short nDim=3;
  double Prandtl_Turb=.9,heat_flux_jac_i[5]{};
  void SetHeatFluxJacobian(const double* val_Mean_PrimVar,const double* val_Mean_SecVar,
      double val_eddy_viscosity,double val_thermal_conductivity,
      double val_heat_capacity_cp,double val_dist_ij)
'''+heat+'\n};\n'
    generated=folder/'NativeGeneralJacobian.hpp'
    generated.write_text(header,encoding='utf8')
    own=['scripts/verifySU2IdealAirAdapter.py','scripts/validateIdealAir.py',
         'tests/testSU2IdealAirAdapter.cpp','third_party/ideal_air/CSU2IdealAir.hpp',
         'third_party/ideal_air/SU2EnergyState.hpp','third_party/ideal_air/IdealAir.hpp',
         'docs/validation/v154/PLAN.md']
    files=[ROOT/name for name in own]+[generated,archive,BASE/'toolchain/w64devkit/bin/g++.exe',
        SOURCE/'SU2_CFD/include/fluid/CFluidModel.hpp',SOURCE/'SU2_CFD/src/fluid/CFluidModel.cpp',
        SOURCE/'SU2_CFD/include/fluid/CSutherland.hpp']
    files += [p for p in BUILD.rglob('*') if p.suffix in ('.obj','.a')]
    locks={p.relative_to(ROOT).as_posix():sha(p) for p in files}
    save(folder/'source_lock.json',locks)
    save(folder/'native_function.json',dict(source='SU2_CFD/src/numerics/CNumerics.cpp',
        signature=signature,bodySHA256=hashlib.sha256(fragment.encode()).hexdigest(),
        heatBodySHA256=hashlib.sha256(heat.encode()).hexdigest(),
        sourceCommit='7048fd35f56d4bf52e586429d6746a36a10255f8',numericBodyUnmodified=True))
    for name in own:
        path=folder/'executed_sources'/name
        path.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(str(ROOT/name),str(path))
    executable=folder/'adapter.exe'
    command=[str(BASE/'toolchain/w64devkit/bin/g++.exe'),'-std=c++17','-O2','-fopenmp',
        '-I'+str(SOURCE/'SU2_CFD/include/fluid'),'-I'+str(BUILD/'Common/src'),
        '-I'+str(SOURCE/'externals/eigen'),'-I'+str(SOURCE/'externals/mel'),'-I'+str(folder),
        str(ROOT/'tests/testSU2IdealAirAdapter.cpp'),'-Wl,--start-group',
        str(BUILD/'SU2_CFD/src/libSU2core.a'),str(BUILD/'Common/src/libSU2Common.a'),
        '-Wl,--end-group','-o',str(executable)]
    try:
        run(folder,'compile',command,120)
        run(folder,'test',[str(executable),str(folder/'results.csv')],60)
        table=np.genfromtxt(str(folder/'results.csv'),delimiter=',',names=True)
        errors={key:float(max(table[key])) for key in table.dtype.names if key.endswith('error')}
        report=dict(passed=True,cases=len(table),maximumErrors=errors,
            rejectedStates=int(sum(table['atomic_rejections'])),executableSHA256=sha(executable),
            bytes=executable.stat().st_size,newQualifiedFlowArea=0)
        save(folder/'report.json',report)
        print(json.dumps(report,indent=2))
    finally:
        stable={name:sha(ROOT/name)==expected for name,expected in locks.items()}
        save(folder/'source_stability.json',stable)
        assert all(stable.values()),'Source or library changed during adapter validation'


if __name__=='__main__':
    main()
