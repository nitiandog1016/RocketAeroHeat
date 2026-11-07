"""唯一已登记真实显式候选；无参数扫描，无自动追加。"""
from pathlib import Path
import json,shutil,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from runPreparedEilmerCase import execute
base=Path(__file__).resolve().parent
dest=base/'source_at_run';dest.mkdir(exist_ok=False)
for name in ['runPreparedEilmerCase.py','runEilmerProcess.py','eilmerExplicitProfile.py','eilmerExplicitAdmission.py',
             'auditEilmerCompletion.py','auditEilmerExplicitCompletion.py','readEilmerNative.py']:
    shutil.copyfile(ROOT/'scripts'/name,dest/name)
r=execute(base/'real01','run',base/'real01_run',660,1024,base/'real01_prepare/report.json')
print(json.dumps({k:r.get(k) for k in ['passed','elapsedSeconds','error','nativeCompletion']}))
raise SystemExit(0 if r['passed'] else 1)
