"""拒绝并列自由分支、错误基态和伪来源被标成温度投影实际结果。"""
import json
from pathlib import Path
from buildRocketTemperatureWorkbench import ROOT,OUT,COMP,bind
from screenRocketFullFLeastSquares import sha
study=COMP/'temperature_projection_comparison.json';candidate=COMP/'projected/direction.mat'
state=ROOT/'outputs/v83_validation/krylov28/wall_velocity_krylov_evidence.mat'
valid=bind(study,candidate,state,'projected');assert valid['constraintApplied']
rejected=[]
for case in ('krylov26','krylov27'):
    try:bind(study,candidate,ROOT/'outputs/v83_validation'/case/'wall_velocity_krylov_evidence.mat','projected')
    except AssertionError:rejected.append(case)
bad=OUT/'bad_temperature_study_source.json';assert not bad.exists();r=json.loads(study.read_text('utf8'));r['sourceBinding']['sourceFileSHA256']='0'*64;bad.write_text(json.dumps(r),encoding='utf8')
try:bind(bad,candidate,state,'projected')
except AssertionError:rejected.append('wrong-source-SHA')
assert len(rejected)==3
out=dict(status='passed-wrong-base-parallel-state-and-source-SHA-rejections',rejected=rejected,
         studySHA256=sha(study),candidateSHA256=sha(candidate),scriptSHA256=sha(__file__),newProductionSolvedAreaM2=0)
target=OUT/'view_binding_rejections.json';assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
