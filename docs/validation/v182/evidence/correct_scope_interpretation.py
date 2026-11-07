"""纠正首版解释的过强断言，保持原数值和首版报告；不重跑原生程序。"""
from pathlib import Path
import json,time
base=Path(__file__).resolve().parent;start=time.monotonic()
original=json.loads((base/'scope01/report.json').read_text())
differences={n:max(v['matrixPostResidual_vs_originalOrder2']) for n,v in original['records'].items()}
assert max(differences.values())>1e-4
r=dict(passed=False,classificationCompleted=True,supersedesInterpretation='scope01/report.json',
       rejectedOriginalStatement='observed post-matrix residual matches original order 2',
       measuredPostMatrixVsOriginalOrder2MaxByCase=differences,
       supportedSourceFact='evaluate_jacobian saves/restores GlobalConfig.interpolation_order rather than block phase order',
       supportedNumericalFact='The huge post-matrix real Jv difference is accounted for by its mismatched R0 offset to about 5.2e-16; this algebra does not isolate every state/cache side effect.',
       validNativeComparison='Before-matrix mixed-order Jv vs independent central difference: late real <=3.94e-6, complex <=1.83e-9; initial both approximately0.2175, so initial gate remains failed.',
       productionCauseOfV181Proven=False,uniqueCauseUnknown=True,
       routeDecision='Close current real-NK full-domain route; select native explicit local pseudo-time marching for next separately budgeted entry.',
       nativeRuns=0,physicalQualification=False,newQualifiedAreaM2=0,elapsedSeconds=time.monotonic()-start)
(base/'scope02.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
