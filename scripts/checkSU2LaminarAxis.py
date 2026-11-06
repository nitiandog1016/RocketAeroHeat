"""模块用途：逐字抽取已构建源项，先检验真实轴线，再决定能否进入层流鼻部计算。"""
from pathlib import Path
import argparse
import csv
import json
import sys

import numpy as np

from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import BIN, run, save

SOURCE = ROOT/'tools/runtime/su2-thermal-build/source'
BUILD = ROOT/'outputs/v159_validation/build01'
RUNTIME = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-axisymmetric-v1'
FILE = 'SU2_CFD/src/numerics/flow/flow_sources.cpp'
HEADER = 'SU2_CFD/include/numerics/flow/flow_sources.hpp'
AUX = 'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl'


def method(text, signature):
    """模块一：仅截取完整函数，不替换算式，不把 Python 重写当原生证明。"""
    begin = text.index(signature)
    opening = text.index('{', begin)
    level = 1
    end = opening+1
    while level:
        level += (text[end] == '{')-(text[end] == '}')
        end += 1
    return text[begin:end]


STUB = r'''
#include <cmath>
#include <cstdio>
#include <algorithm>
using su2double = double;
constexpr double EPS = 1e-16, TWO3 = 2.0/3.0;
struct CConfig {};
struct CNumerics {
  template<class T=double> struct ResidualType {
    const double* residual;
    ResidualType(const double* r, double**, void*) : residual(r) {}
  };
};
// 最小存储壳，仅提供逐字原函数读取的成员；不是整个 SU2 类/离散求解器。
struct CSourceAxisymmetric_Flow : CNumerics {
  unsigned short nVar=4, nDim=2;
  bool implicit=true, viscous=true, rans=false;
  double yinv=0.0, Volume=.7, turb_ke_i=0, Prandtl_Turb=.9;
  double Coord_i[2]{}, U_i[4]{}, V_i[12]{}, V_j[12]{}, S_i[2]{};
  double PrimVar_Grad_i[3][2]{}, AuxVar_Grad_i[3][2]{};
  double residual[4]{}, matrix[4][4]{};
  double* jacobian[4]{matrix[0], matrix[1], matrix[2], matrix[3]};
  void ResidualDiffusion();
};
struct CSourceGeneralAxisymmetric_Flow : CSourceAxisymmetric_Flow {
  ResidualType<> ComputeResidual(const CConfig*);
};
'''

DRIVER = r'''
void fill(CSourceGeneralAxisymmetric_Flow& x, double r, int mode) {
  // 光滑偶/奇场：u=10+a*r²，v=c*r，T=300+b*r²；轴线极限有解析值。
  const double a=mode==2 ? 2.0 : 0.0;
  const double b=mode==2 ? 5.0 : 0.0;
  const double c=mode==1 ? .3 : 0.0;
  const double rho=1.2, u=10+a*r*r, v=c*r, mu=2e-5, k=.03;
  const double p=100000, e=250000;
  x.Coord_i[1]=r;
  x.U_i[0]=rho; x.U_i[1]=rho*u; x.U_i[2]=rho*v;
  x.U_i[3]=rho*(e+.5*(u*u+v*v));
  x.V_j[3]=p; x.V_i[7]=mu; x.V_i[8]=0; x.V_i[9]=k; x.V_i[10]=1005;
  x.S_i[0]=p/rho; x.S_i[1]=.4*rho;
  x.PrimVar_Grad_i[0][1]=2*b*r;
  x.PrimVar_Grad_i[1][1]=2*a*r;
  x.PrimVar_Grad_i[2][1]=c;
  // 原辅助变量为 mu*v/r、mu*v²/r、mu*u*v/r。第二项径向导数不能遗漏。
  for(auto& row : x.AuxVar_Grad_i) for(auto& v : row) v=0;
  x.AuxVar_Grad_i[1][1]=mu*c*c;
}
int main() {
  std::puts("mode,viscous,preceding_radius,radius,R0,R1,R2,R3");
  for(int mode : {1,2}) for(bool viscous : {false,true}) {
    for(double before : {0.0,.1,.2}) {
      CSourceGeneralAxisymmetric_Flow x; x.viscous=viscous;
      if(before>0) {fill(x,before,mode); x.ComputeResidual(nullptr);}
      fill(x,0,mode); x.ComputeResidual(nullptr);
      std::printf("%d,%d,%.17g,0,%.17g,%.17g,%.17g,%.17g\n",mode,viscous,before,
                  x.residual[0],x.residual[1],x.residual[2],x.residual[3]);
    }
    CSourceGeneralAxisymmetric_Flow x; x.viscous=viscous;
    for(double radius : {.01,.001,.0001}) {
      fill(x,radius,mode); x.ComputeResidual(nullptr);
      std::printf("%d,%d,0,%.17g,%.17g,%.17g,%.17g,%.17g\n",mode,viscous,radius,
                  x.residual[0],x.residual[1],x.residual[2],x.residual[3]);
    }
  }
}
'''


def check(folder):
    """模块二：绑定准确源树/可执行文件，编小检查并保留通过与失败，绝不启动 CFD。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    folder.mkdir(parents=True, exist_ok=False)
    (folder/'executed_check.py').write_bytes(Path(__file__).read_bytes())
    runtime = read(RUNTIME/'runtime.json')
    tree_file = BUILD/'source_tree.json'
    assert sha(tree_file) == runtime['sourceTreeSHA256']
    tree = read(tree_file)
    paths = [SOURCE/FILE, SOURCE/HEADER, SOURCE/AUX]
    locks = {str(p): sha(p) for p in paths+[tree_file, RUNTIME/'runtime.json', Path(__file__).resolve(),
                                           ROOT/'docs/validation/v169/PLAN.md']}
    for p in paths:
        assert sha(p) == tree[str(p)], '当前共享源不是该运行包的准确构建源'
        (folder/p.name).write_bytes(p.read_bytes())
    for entry in runtime['files']:
        assert sha(RUNTIME/entry['path']) == entry['sha256']
    exe = RUNTIME/'bin/SU2_CFD.exe'
    locks[str(exe)] = sha(exe)
    source = paths[0].read_text(encoding='utf8')
    diffusion = method(source, 'void CSourceAxisymmetric_Flow::ResidualDiffusion()')
    general = method(source, 'CNumerics::ResidualType<> CSourceGeneralAxisymmetric_Flow::ComputeResidual(')
    (folder/'axis_methods.cpp').write_text(STUB+'\n'+diffusion+'\n'+general+'\n'+DRIVER, encoding='utf8')
    run(folder, 'compile', [str(BIN/'g++.exe'), '-std=c++17', '-O2', str(folder/'axis_methods.cpp'),
                           '-o', str(folder/'axis_methods.exe')], 30)
    run(folder, 'evaluate', [str(folder/'axis_methods.exe')], 15)
    rows = np.genfromtxt(str(folder/'evaluate.log'), delimiter=',', names=True)
    result = []
    # 独立极限：v=c*r 的径向源应为零；u/T 偶二次项产生非零黏性/导热极限。
    for mode in (1, 2):
        for viscous in (False, True):
            selected = rows[(rows['mode'] == mode) & (rows['viscous'] == int(viscous))]
            actual = np.column_stack([selected['R%d'%j] for j in range(4)])
            expected = []
            for radius in selected['radius']:
                a, b, c = (0.,0.,.3) if mode == 1 else (2.,5.,0.)
                rho, mu, k, vol = 1.2, 2e-5, .03, .7
                u, v = 10+a*radius**2, c*radius
                h = 250000+.5*(u*u+v*v)+100000/rho
                q = vol*np.array([rho*c, rho*u*c, rho*c*c*radius, rho*h*c])
                if viscous:
                    q -= vol*np.array([0.,2*mu*a,0.,2*mu*a*u+2*k*b])
                expected.append(q)
            expected = np.array(expected)
            axis = selected['radius'] == 0
            error = np.max(np.abs(actual-expected)/np.maximum(1., abs(expected)), axis=1)
            order_error = float(np.max(np.ptp(actual[axis], axis=0)))
            result.append(dict(mode=mode, viscous=viscous, axisMaximumScaledError=float(max(error[axis])),
                positiveRadiusMaximumScaledError=float(max(error[~axis])),
                axisCallOrderAbsoluteDifference=order_error,
                axisLimitPassed=bool(max(error[axis]) < 1e-12), callOrderPassed=order_error < 1e-12,
                positiveRadiusControlPassed=bool(max(error[~axis]) < 1e-12),
                actualAxis=actual[axis].tolist(), expectedAxis=expected[axis].tolist()))
    stable = {p: sha(Path(p)) == h for p,h in locks.items()}
    assert all(stable.values())
    controls = all(r['positiveRadiusControlPassed'] for r in result) and all(
        r['axisLimitPassed'] for r in result if not r['viscous'])
    passed = all(r['axisLimitPassed'] and r['callOrderPassed'] for r in result)
    save(folder/'report.json', dict(sourceBound=True, harness='verbatim-methods-minimal-storage-not-full-class',
        controlsPassed=controls, axisymmetricLaminarEntryPassed=passed, results=result,
        inputSHA256=locks, inputsStable=stable, newCFDRuns=0, newQualifiedArea_m2=0.,
        decision='eligible-for-small-native-check' if passed else 'reject-current-axisymmetric-nose-entry',
        limitations='local source check, not whole-discretization proof or cause of old positive-radius/3D failures'))
    assert controls, '控制组失败，不能据此判定引擎缺项'
    print(json.dumps(dict(controlsPassed=controls, axisymmetricLaminarEntryPassed=passed)))
    return 0 if passed else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    args = parser.parse_args()
    raise SystemExit(check(args.folder))
