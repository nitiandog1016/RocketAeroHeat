from pathlib import Path
import hashlib,json,shutil
root=Path.cwd();base=root/'outputs/v168_validation';dest=root/'docs/validation/v168/evidence'
assert not dest.exists()
files=list(base.glob('*_process.json'))+list(base.glob('*.log'))
for name in ['review01','review02','benchmark01','verify01','verify02','verify03']:
    files += [p for p in (base/name).rglob('*') if p.is_file() and p.suffix in ('.py','.md','.json','.txt','.csv','.log')]
files += [base/'sources/NACA_TN_3344_source.json']
files += [Path(__file__).resolve()]
records=[]
for p in sorted(set(files)):
    p.resolve().relative_to(root.resolve());r=p.relative_to(root);target=dest/r
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(target))
    sha=hashlib.sha256(p.read_bytes()).hexdigest();assert hashlib.sha256(target.read_bytes()).hexdigest()==sha
    records.append(dict(source=r.as_posix(),path=target.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=sha))
fig=root/'docs/validation/v168/figures';fig.mkdir(exist_ok=False)
for name in ['wall_heat_review.png','naca3344_digitized.png']:
    p=base/'verify03'/name;target=fig/name;shutil.copyfile(str(p),str(target))
    records.append(dict(source=p.relative_to(root).as_posix(),path=target.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
result=dict(files=records,fileCount=len(records),bytes=sum(x['bytes'] for x in records),uniqueDataDeleted=0,fullProjectPackages=0,largeCFDOutputsCopied=0,originalPDFPublished=False)
(root/'docs/validation/v168/evidence_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in result.items() if k!='files'}))
