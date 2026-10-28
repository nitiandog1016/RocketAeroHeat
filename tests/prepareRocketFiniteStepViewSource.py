"""为既有导出器建立真实接受检查点的只读展示别名，不新增计算或接受历史。"""
from pathlib import Path
import json,os
from recordRocketFiniteStepProgress import sha
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v84_validation'
full=json.loads((OUT/'full02/full_comparison.json').read_text('utf8'))
candidates=[(step['main']['rmsMerit'],step['main']['maximumResidual'],branch['method'],step) for branch in full['branches'] for step in branch['steps'] if step['main'] and step['main']['acceptedNonlinearIterations']>0]
best=min(candidates,key=lambda x:(x[0],x[1]));step=best[3]
# 只在所有候选均无更低RMS且更低峰值的反例时保留明确选择；记录实际对照身份。
source=ROOT/step['resultFolder'];source.resolve().relative_to(ROOT.resolve())
alias=ROOT/'outputs/v83_validation/krylov25';assert not alias.exists();alias.mkdir()
for name in ['wall_velocity_krylov_evidence.mat','wall_velocity_krylov_summary.json']:
 os.link(source/name,alias/name);assert sha(source/name)==sha(alias/name)
binding={'status':'byte-identical-display-alias-of-actual-accepted-state','sourceFolder':str(source.relative_to(ROOT)),
 'displayCase':'krylov25','selectedMethod':best[2],'selectedRMS':best[0],'selectedMaximum':best[1],
 'fileSHA256':sha(alias/'wall_velocity_krylov_evidence.mat'),'newFlowAcceptances':0,'newProductionSolvedAreaM2':0}
(OUT/'view_source_binding.json').write_text(json.dumps(binding,indent=2)+'\n',encoding='utf8');print(json.dumps(binding))
