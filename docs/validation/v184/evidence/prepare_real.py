"""唯一真实显式候选准备，复用审核网格与原均匀初态。"""
from pathlib import Path
import json,shutil,sys
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from prepareEilmerDevelopmentCase import prepare
base=Path(__file__).resolve().parent;source=ROOT/'outputs/v178_validation/real01'
dest=base/'source_at_prepare';dest.mkdir(exist_ok=False)
for name in ['prepareEilmerDevelopmentCase.py','prepareEilmerRocketCase.py','eilmerExplicitProfile.py','eilmerExplicitAdmission.py','readEilmerNative.py']:
    shutil.copyfile(ROOT/'scripts'/name,dest/name)
r=prepare(source,ROOT/'outputs/v178_validation/real01_prepare/report.json',ROOT/'outputs/v178_validation/real01_geometry.json',
          base/'real01',base/'real01_prepare',json.loads((base/'condition.json').read_text()),120,1024)
print(json.dumps({k:r.get(k) for k in ['passed','elapsedSeconds','error','nativeConfigChangedKeys']}))
raise SystemExit(0 if r['passed'] else 1)
