"""唯一完整域新存储任务；每阶段保存准确执行源，预算与候选登记一致。"""
from pathlib import Path
import json
import shutil
import sys

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from prepareEilmerDevelopmentCase import prepare
from runPreparedEilmerCase import execute

base=Path(__file__).resolve().parent;phase=sys.argv[1]
snapshot=base/('source_at_'+phase);snapshot.mkdir(exist_ok=False)
for name in ['prepareEilmerDevelopmentCase.py','prepareEilmerRocketCase.py','eilmerExplicitProfile.py',
             'eilmerExplicitAdmission.py','eilmerDenseRuntime.py','runPreparedEilmerCase.py','runEilmerProcess.py',
             'readEilmerNative.py','auditEilmerExplicitCompletion.py','auditEilmerCompletion.py','inspectEilmerRun.py']:
    shutil.copyfile(ROOT/'scripts'/name,snapshot/name)
if phase=='prepare':
    result=prepare(ROOT/'outputs/v178_validation/real01',ROOT/'outputs/v178_validation/real01_prepare/report.json',
                   ROOT/'outputs/v178_validation/real01_geometry.json',base/'real01',base/'real01_prepare',
                   json.loads((base/'condition.json').read_text()),120,1024)
else:
    assert phase=='run'
    result=execute(base/'real01','run',base/'real01_run',660,1024,base/'real01_prepare/report.json')
print(json.dumps({k:result.get(k) for k in ['passed','elapsedSeconds','error']}))
raise SystemExit(0 if result['passed'] else 1)
