"""唯一六步显式小例；固定预算，不做失败重算或参数扫描。"""
from pathlib import Path
import json,shutil,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from runPreparedEilmerCase import execute
BASE=Path(__file__).resolve().parent
dest=BASE/'source_at_run';dest.mkdir(exist_ok=False)
for name in ['runPreparedEilmerCase.py','auditEilmerCompletion.py','auditEilmerExplicitCompletion.py',
             'eilmerExplicitProfile.py','readEilmerNative.py','runEilmerProcess.py']:
    shutil.copyfile(ROOT/'scripts'/name,dest/name)
r=execute(BASE/'small01','run',BASE/'small01_run',30,128,BASE/'small01_prepare/report.json')
print(json.dumps({k:r.get(k) for k in ['passed','elapsedSeconds','error','nativeCompletion']}))
raise SystemExit(0 if r['passed'] else 1)
