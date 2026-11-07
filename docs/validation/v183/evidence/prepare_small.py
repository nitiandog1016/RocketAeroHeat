"""只准备一次固定显式小例；准确源和原生失败均保留。"""
from pathlib import Path
import json,shutil,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from prepareEilmerDevelopmentCase import prepare
BASE=Path(__file__).resolve().parent
source=ROOT/'outputs/v178_validation/native03'
c=json.loads((source/'case_manifest.json').read_text())['condition']
c.pop('maxNewtonSteps');c['explicit']=dict(profile='rah-explicit-local-v1',maxSteps=6,checkpointStep=3)
(BASE/'condition.json').write_text(json.dumps(c,indent=2)+'\n')
dest=BASE/'source_at_prepare';dest.mkdir(exist_ok=False)
for name in ['prepareEilmerDevelopmentCase.py','prepareEilmerRocketCase.py','eilmerExplicitProfile.py',
             'runPreparedEilmerCase.py','auditEilmerExplicitCompletion.py','readEilmerNative.py','auditEilmerCompletion.py']:
    shutil.copyfile(ROOT/'scripts'/name,dest/name)
r=prepare(source,ROOT/'outputs/v178_validation/native03_prepare/report.json',
          ROOT/'outputs/v178_validation/native03_geometry.json',BASE/'small01',BASE/'small01_prepare',c,60,128)
print(json.dumps({k:r.get(k) for k in ['passed','elapsedSeconds','error','nativeConfigChangedKeys']}))
raise SystemExit(0 if r['passed'] else 1)
