"""复核已保存原生量与准确源码：区分插值阶副作用、失配基准列和有效同阶导数。"""
from pathlib import Path
import hashlib,json,re,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];base=ROOT/'outputs/v182_validation';start=time.monotonic()
out=base/'scope01';out.mkdir(exist_ok=False)
records={}
for name in ['real-state0','complex-state0','real-state2','complex-state2']:
    a=np.genfromtxt(base/name/'operator.csv',delimiter=',',names=True)
    log=(base/(name+'-process')/'stdout.txt').read_text()
    sigma=float(re.search(r'RAH_SIGMA (\S+)',log)[1])
    def err(x,y):return [float(np.linalg.norm((x-y)[a['variable']==i])/max(1.,np.linalg.norm(y[a['variable']==i]))) for i in range(7)]
    d=dict(matrixPostResidual_vs_originalOrder2=err(a['R_after_matrix'],a['R_high_before']),
           validMixedJv_vs_central=err(a['J_mixed'],a['J_central']),
           sparseJv_vs_central=err(a['J_sparse'],a['J_central']))
    if name.startswith('real'):
        # 只做代数归因核对，不冒充已经运行的修正Jv。
        accounted=a['J_mixed']+(a['R']-a['R_after_matrix'])/sigma
        d['mismatchedR0ExplainsPostMatrixJv']=err(a['J_native'],accounted)
    records[name]=d
src=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/fluidblock.d');lines=src.read_text().splitlines()
snippet='\n'.join(f'{i+1}: {lines[i]}' for i in range(1478,1543))+'\n';(out/'fluidblock-order-scope.txt').write_text(snippet)
r=dict(passed=True,scope='classification-of-existing-diagnostic-columns',records=records,
       sourceSHA256=hashlib.sha256(src.read_bytes()).hexdigest(),
       observedOrderRestoration='evaluate_jacobian saves GlobalConfig.interpolation_order, not the current block phase order; observed post-matrix residual matches original order 2.',
       invalidColumn='Post-matrix J_native in real diagnostics used order-2 R0 after the matrix reset while Jv evaluates order 1. Do not call its enormous error a physical derivative error.',
       stillValidColumns=['J_mixed before matrix; native production mixed-order R0','J_sparse at explicit Jacobian order','J_central at matching order','R/U restoration observations'],
       correctedNativeRerun=False,productionCauseOfV181Proven=False,initialStateDerivativeGateStillFailed=True,
       realNKRouteStatus='closed-for-full-domain; no proven minimal production fix closes all required checks',
       newNativeRuns=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
