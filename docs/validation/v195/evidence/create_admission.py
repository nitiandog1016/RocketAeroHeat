"""检查已通过的真实初场接口/监督器后，绑定唯一400步候选。"""
from pathlib import Path
import json
import sys
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'scripts'))
from hybridRocketEuler import BUDGET,sha,save,check_admission
base=root/'outputs/v195_validation'
for path in [base/'interface01/report.json',base/'guard01/report.json']:
    assert json.loads(path.read_text())['passed']
paths=[base/'interface01/report.json',base/'guard01/report.json',base/'preparation01/report.json',
       root/'docs/validation/v195/PLAN.md']
paths += [root/'scripts'/name for name in ['hybridRocketEuler.py','runEilmerProcess.py','readEilmerNative.py',
           'prepareHybridRocketEulerNative.py','eilmerAirThermo.py','hybridExternalState.py','hybridEdgeSampling.py']]
runtime=root/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json'
a=dict(schema='RAH.HybridRocketEulerAdmission.v1',profile='full-wall-euler-400-v1',budget=BUDGET,
       temperatureRangeK=[200.,6000.],case=str(base/'rocket01'),preparation=str(base/'preparation01/report.json'),
       record=str(base/'run01'),runtimeSHA256=sha(runtime),evidenceSHA256={str(p):sha(p) for p in paths})
assert not (base/'admission.json').exists();check_admission(a);save(base/'admission.json',a)
print('Unique 400-step admission bound; no execution in this script.')
