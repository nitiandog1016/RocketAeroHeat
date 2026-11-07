"""有限核对已执行包的准确源；归档时间尺度/退化信号/保存循环，不称已证明唯一根因。"""
from pathlib import Path
import hashlib,json,time
start=time.monotonic();root=Path('/mnt/e/RocketAeroHeat');out=root/'outputs/v195_validation/source_review01'
out.mkdir(exist_ok=False)
build=json.loads((root/'outputs/v187_validation/build01/report.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
src=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source/src')
spec=[(src/'lmr/fluidfvcell.d',[(463,490),(1000,1083)]),
      (src/'geom/elements/properties.d',[(326,345)]),
      (src/'lmr/turbulence.d',[(109,112)]),
      (src/'lmr/timemarching.d',[(80,94),(603,615)]),
      (root/'tools/runtime/eilmer-explicit-workspace-v1/overlay/src/lmr/fluidblock.d',[(1200,1255)])]
rows=[]
for p,ranges in spec:
    digest=sha(p);assert build['inputSHA256'][str(p)]==digest
    lines=p.read_text().splitlines()
    rows.append(dict(path=str(p),sha256=digest,boundToActualBuild=True,
        excerpts=[dict(firstLine=a,lastLine=b,text='\n'.join('%d:%s'%(i,lines[i-1]) for i in range(a,b+1))) for a,b in ranges]))
result=dict(passed=True,sourceCommit=build['sourceCommit'],actualBuildReceiptSHA256=sha(root/'outputs/v187_validation/build01/report.json'),
    sources=rows,findings=['geometry source assigns all tetrahedral lengths from computed L_min',
        'signal_frequency reads selected normals, sound speed and directional lengths; no actual runtime values observed',
        'none turbulence frequency is 1e-99; 0.3 divided by that value matches observed 3e98, but this is not proof of why convective signal disappeared',
        'local dt is capped separately; global dt_allow is returned from unbounded cfl/signal',
        'snapshot scheduler adds dt from zero in a while loop; huge now makes this impractical'],
    uniqueRuntimeCauseEstablished=False,newBuilds=0,newCFDRuns=0,elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
(out/'executed_review.py').write_bytes(Path(__file__).read_bytes());print(json.dumps(dict(passed=True,exactSources=len(rows),elapsedSeconds=result['elapsedSeconds'])))
