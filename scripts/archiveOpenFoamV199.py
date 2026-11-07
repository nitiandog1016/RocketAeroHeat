"""v199公开小证据清单：原字节绑定，体网格/体场/运行库保留本地。"""
from pathlib import Path
import hashlib,json,shutil,time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v199_validation';DEST=ROOT/'docs/validation/v199/evidence'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False)
    # 模块一：只复制明确的小证据类型。没有扩展名的体场、网格及npz不在集合内。
    suffixes={'.json','.txt','.log','.py','.C','.H','.md','.csv','.png'}
    special={'controlDict','fvSchemes','fvSolution','physicalProperties'}
    rows=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        if p.suffix not in suffixes and p.name not in special:continue
        # 合成100000行日志只留摘要/准确测试源，不上传重复写入文本。
        if p.relative_to(BASE).as_posix()=='native_log01/native/stdout.txt':continue
        q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==sha(p)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),
            bytes=p.stat().st_size,sha256=sha(p)))
    # 模块二：本版关键最终源码与目标同证据绑定；准确执行源另有原运行快照。
    final=['AGENTS.md','docs/PROJECT_GOALS.md','docs/HYBRID_SYNCHRONOUS_CONTROL_V199.md',
        'docs/validation/v199/PLAN.md','docs/validation/v199/DEVELOPMENT_PLAN.md','docs/validation/v199/NEXT_DECISION.md',
        'scripts/runOpenFoamProcess.py','scripts/openFoamConservation.py','scripts/auditOpenFoamDevelopment.py',
        'scripts/openFoamSynchronousControl.py','scripts/reviewOpenFoamSynchronous.py','scripts/auditOpenFoamSynchronous.py',
        'scripts/archiveOpenFoamV199.py','tests/testOpenFoamConservation.py','tests/testOpenFoamObservation.py',
        'tests/testOpenFoamNativeLog.py']
    for name in final:
        p=ROOT/name;q=DEST/'final_source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    report=dict(schema='RAH.V199Evidence.v1',files=rows,count=len(rows),bytes=sum(x['bytes'] for x in rows),
        originalFilesUnmodified=True,protectedFiles=29,largeAssetsUploaded=False,uniqueDeleted=0,
        outputLogicalBytes=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file()),
        elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v199/evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['count','bytes','outputLogicalBytes','protectedFiles','elapsedSeconds']}))


if __name__=='__main__':main()
