"""只读确认实数研究结果携带具体数值方案及失败资格，避免和复数方法混称。"""
from pathlib import Path
import csv,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import export_loads
base=ROOT/'outputs/v181_validation';out=base/'identity01';out.mkdir(exist_ok=False)
start=time.monotonic();case=base/'small01'
r=export_loads(case,case/'lmrsim/loads/0003/blk-0000-bndry-0-wall.dat',out/'wall.csv',base/'small01_run/report.json')
rows=list(csv.DictReader((out/'wall.csv').open()))
assert len(rows)==7 and all('rah-nk-real-v1' in x['method'] and 'comparison failed' in x['method'] and x['physicsQualified']=='False' for x in rows)
assert r['compatibilityStatus']=='small-case-comparison-failed' and r['newQualifiedAreaM2']==0
r.update(scope='Export identity only; original comparison remains failed',newNativeRuns=0,
         sourceSHA256=hashlib.sha256((ROOT/'scripts/readEilmerNative.py').read_bytes()).hexdigest(),elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
