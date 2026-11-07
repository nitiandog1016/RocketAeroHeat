from pathlib import Path
import hashlib,json,shutil,time
root=Path.cwd();base=root/'outputs/v178_validation';dest=root/'docs/validation/v178/evidence'
assert not dest.exists();dest.mkdir(parents=True)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
protected=json.loads((root/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
checks={p:sha(root/p)==v for p,v in protected.items()};assert all(checks.values())
(base/'protected_final.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf8')
files=[];excluded=[]
for p in sorted(base.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts or p.suffix=='.pyc':continue
    if p.stat().st_size>1024*1024:
        excluded.append(dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=sha(p)));continue
    target=dest/'runs'/p.relative_to(base);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(target))
    files.append(dict(path=target.relative_to(root).as_posix(),source=p.relative_to(root).as_posix(),sha256=sha(target),bytes=target.stat().st_size))
for name in ['reference01','reference_oracle01']:
    for p in sorted((root/'outputs/v177_validation'/name).rglob('*')):
        if not p.is_file():continue
        t=dest/'reference'/name/p.relative_to(root/'outputs/v177_validation'/name);t.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(t))
        files.append(dict(path=t.relative_to(root).as_posix(),source=p.relative_to(root).as_posix(),sha256=sha(t),bytes=t.stat().st_size))
sources=['runEilmerCase.m','tests/testEilmerMatlabEntry.m','scripts/prepareEilmerRocketCase.py','scripts/readEilmerNative.py','scripts/runEilmerProcess.py','scripts/runPreparedEilmerCase.py','tests/testEilmerRocketPreparation.py','tests/testEilmerNativeInterface.py','tests/testEilmerProcessSupervision.py','tests/verifyEilmerWallGeometry.py','scripts/prepareEilmerCylinderFlareReference.py','tests/verifyEilmerCylinderFlareReference.py']
for name in sources:
    p=root/name;t=dest/'source_versions'/(p.name+'.'+sha(p)+'.raw');t.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(p),str(t))
    files.append(dict(path=t.relative_to(root).as_posix(),source=name,sha256=sha(t),bytes=t.stat().st_size))
report=dict(files=files,filesCount=len(files),totalSmallBytes=sum(x['bytes'] for x in files),largeLocalOnly=excluded,protectedFiles=len(checks),allProtectedStable=True,uniqueFilesDeleted=0,newQualifiedAreaM2=0)
(root/'docs/validation/v178/evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in report.items() if k not in ('files','largeLocalOnly')}))
print(json.dumps(excluded))
