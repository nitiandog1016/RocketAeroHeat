"""保存显式原生输出/残差的准确源码证据，不运行或修改引擎。"""
from pathlib import Path
import hashlib,json,time
start=time.monotonic();out=Path(__file__).resolve().parent/'source_contract';out.mkdir(exist_ok=False)
root=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr')
ranges={'lua-modules/output.lua':[(8,36),(70,135)],'timemarching.d':[(115,123),(166,175),(649,658),(731,751),(789,807),(880,897)],
        'fluidblock.d':[(1004,1048)],'commands/runsim.d':[(342,367)]}
files=[]
for rel,spans in ranges.items():
    source=root/rel;raw=source.read_bytes();lines=raw.decode('utf8').splitlines()
    text='\n\n'.join('\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1)) for a,b in spans)+'\n'
    dest=out/(rel.replace('/','__')+'.txt');dest.write_text(text)
    files.append(dict(source=str(source),sha256=hashlib.sha256(raw).hexdigest(),spans=spans,excerpt=dest.name))
r=dict(sourceCommit='f53f4609a0331d48efee69a4e4f3c3598378cc03',files=files,newNativeRuns=0,
       notes=['Direct lmr-run uses default verbosity; wrapper introduction is conditional, not mandatory.',
              'initLoadsFiles initializes only; checkpoint is first native load index 0.',
              'Transient Linf uses dUdt[0], not delta-U; RK3 stage-0 residual timing differs from final saved state.',
              'Printed mixed-variable L2 is unscaled and not the NK relative-residual metric; no physical convergence claimed.'],
       elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(files=len(files),elapsedSeconds=r['elapsedSeconds'])))
