"""唯一新存储包小例：复用原几何/初态，每个阶段各有父监督和准确执行源。"""
from pathlib import Path
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from runEilmerProcess import sha
from prepareEilmerDevelopmentCase import prepare
from runPreparedEilmerCase import execute

base = Path(__file__).resolve().parent
phase = sys.argv[1]
snapshot = base / ('source_at_' + phase)
snapshot.mkdir(exist_ok=False)
for name in ['prepareEilmerDevelopmentCase.py', 'prepareEilmerRocketCase.py', 'eilmerExplicitProfile.py',
             'eilmerDenseRuntime.py', 'runPreparedEilmerCase.py', 'runEilmerProcess.py',
             'readEilmerNative.py', 'auditEilmerExplicitCompletion.py', 'auditEilmerCompletion.py', 'inspectEilmerRun.py']:
    shutil.copyfile(ROOT / 'scripts' / name, snapshot / name)
if phase == 'prepare':
    source = ROOT / 'outputs/v178_validation/native03'
    condition = json.loads((ROOT / 'outputs/v183_validation/condition.json').read_text())
    package = ROOT / 'tools/runtime/eilmer-explicit-dense-v1/runtime.json'
    condition['explicit'].update(profile='rah-explicit-dense-v1', runtimePackage=str(package), runtimePackageSHA256=sha(package))
    (base / 'condition.json').write_text(json.dumps(condition, indent=2)+'\n')
    result = prepare(source, ROOT / 'outputs/v178_validation/native03_prepare/report.json',
                     ROOT / 'outputs/v178_validation/native03_geometry.json', base / 'small01',
                     base / 'small01_prepare', condition, 30, 64)
else:
    assert phase == 'run'
    result = execute(base / 'small01', 'run', base / 'small01_run', 30, 64, base / 'small01_prepare/report.json')
print(json.dumps({k: result.get(k) for k in ['passed', 'elapsedSeconds', 'error']}))
raise SystemExit(0 if result['passed'] else 1)
