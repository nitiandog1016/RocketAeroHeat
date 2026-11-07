"""只把既有LTO对象转为可检查的隔离对象；原对象/安装绝不修改。"""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import run,PREFIX,sha
base=ROOT/'outputs/v182_validation';out=base/'build01';out.mkdir(exist_ok=False)
start=time.monotonic();records=[]
for name in ['lmr-run','lmrZ-run']:
    source=PREFIX/'source/src/lmr'/(name+'.o');target=out/(name+'-native.o')
    report=run(['ld','-r','--plugin='+str(PREFIX/'ldc2-1.41.0-linux-x86_64/lib/LLVMgold-ldc.so'),
                '-plugin-opt=O2',str(source),'-o',str(target)],out,base/(name+'-convert'),min(60,180-(time.monotonic()-start)),256,[source,Path(__file__)])
    records.append(report)
    if not report['passed']:break
    assert target.read_bytes()[:4]==b'\x7fELF','转换仍为bitcode，停止当前复用路线'
r=dict(passed=len(records)==2 and all(x['passed'] for x in records),processes=records,
       elapsedSeconds=time.monotonic()-start,sourceKernelCompiles=0)
(out/'reuse_report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:r[k] for k in ['passed','elapsedSeconds','sourceKernelCompiles']}))
