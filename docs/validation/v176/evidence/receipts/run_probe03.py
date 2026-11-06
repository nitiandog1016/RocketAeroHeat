from pathlib import Path
import sys,json,hashlib,shutil
sys.path.insert(0,str(Path.cwd()/'scripts'))
from validateIdealAir import run
out=Path('outputs/v176_validation')
paths=[r['path'] for r in json.loads((out/'sources02.json').read_text(encoding='utf8'))]+['tests/testAerothermalEvidenceProbe.m']
snapshot=out/'sources03'
assert not snapshot.exists()
manifest=[]
for name in paths:
    source=Path(name);target=snapshot/source
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(source),str(target))
    raw=source.read_bytes();manifest.append(dict(path=name,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)))
(out/'sources03.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')
run(out,'node02',['node','tests/testAerothermalEvidenceView.js'],10)
run(out,'probe03_supervision',[sys.executable,'scripts/runMatlabResearch.py','probe03',"addpath(pwd); addpath('tests'); testAerothermalEvidenceProbe('outputs/v176_validation/probe03',90);",'--folder','outputs/v176_validation','--input','outputs/v104_validation/baseline02/current_case.mat','--input','outputs/v167_validation/main01/result.mat'],170)
