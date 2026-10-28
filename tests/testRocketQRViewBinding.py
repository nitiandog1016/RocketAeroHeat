"""拒绝把错误基态或伪造来源的QR方法贴到当前原面结果。"""
import json
from pathlib import Path
from buildRocketQRWorkbench import ROOT,OUT,bind
from screenRocketFullFLeastSquares import sha

study=OUT/'qr_real01/augmented_qr_study.json';candidate=OUT/'qr_real01/augmented_qr_direction.mat'
state=ROOT/'outputs/v83_validation/krylov26/wall_velocity_krylov_evidence.mat'
valid=bind(study,candidate,state);assert valid['usedTotalDirectionExactlyVerified']
rejected=[]
try:bind(study,candidate,ROOT/'outputs/v83_validation/krylov25/wall_velocity_krylov_evidence.mat')
except AssertionError:rejected.append('wrong-actual-result-state')
bad=OUT/'bad_qr_study_source.json';assert not bad.exists();r=json.loads(study.read_text('utf8'));r['sourceBinding']['sourceFileSHA256']='0'*64
bad.write_text(json.dumps(r),encoding='utf8')
try:bind(bad,candidate,state)
except AssertionError:rejected.append('wrong-source-state-SHA')
assert len(rejected)==2
result=dict(status='passed-wrong-state-and-source-SHA-provenance-rejections',rejected=rejected,
            originalStudySHA256=sha(study),candidateSHA256=sha(candidate),scriptSHA256=sha(__file__),newProductionSolvedAreaM2=0)
target=OUT/'qr_view_binding_rejections.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');print(json.dumps(result))
