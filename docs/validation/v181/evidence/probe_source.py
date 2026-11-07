"""定向读取固定引擎Jv/预条件/实数入口，不执行或修改求解器。"""
from pathlib import Path
import hashlib,json,re,time
start=time.monotonic()
source=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr')
out=Path('/mnt/e/RocketAeroHeat/outputs/v181_validation/source01');out.mkdir(exist_ok=False)
hits=[]
for p in source.glob('*.d'):
    lines=p.read_text().splitlines()
    matches=[i for i,line in enumerate(lines) if re.search(r'preconditionerPerturbation|void evalRealMatVecProd|void computeScaling|useScaling|void scaleConserved|int main\(',line)]
    if not matches:continue
    (out/p.name).write_bytes(p.read_bytes())
    for i in matches:
        hits.append(dict(file=p.name,line=i+1,text=lines[i]))
result=dict(passed=True,hits=hits,sourceSHA256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*.d')},
            initialSearchFailure='WSL rg is unavailable; no solver task was launched',nativeRuns=0,elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
