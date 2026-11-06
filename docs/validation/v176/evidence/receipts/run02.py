from pathlib import Path
import sys, json, hashlib, shutil
sys.path.insert(0,str(Path.cwd()/'scripts'))
from validateIdealAir import run
out=Path('outputs/v176_validation')
original=json.loads((out/'sources01.json').read_text(encoding='utf-8-sig'))
snapshot=out/'sources02'
assert not snapshot.exists()
manifest=[]
for item in original:
    source=Path(item['path'])
    target=snapshot/source
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(str(source),str(target))
    raw=source.read_bytes()
    manifest.append(dict(path=source.as_posix(),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)))
(out/'sources02.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')
run(out,'main02_supervision',[sys.executable,'scripts/runMatlabResearch.py','main02',"addpath(pwd); addpath('tests'); testAerothermalEvidenceWorkbench('outputs/v176_validation/main02',150);",'--folder','outputs/v176_validation','--input','outputs/v104_validation/baseline02/current_case.mat','--input','outputs/v167_validation/main01/result.mat'],300)
