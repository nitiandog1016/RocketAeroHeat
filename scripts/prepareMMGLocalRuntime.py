"""下载固定MMG源，保存版本、许可证与文件身份；不覆盖现有运行时。"""
from pathlib import Path
import json,tarfile,time,urllib.request
from prepareCGALRemeshingRuntime import ROOT,sha
DEST=ROOT/'tools/runtime/mmg/5.8.0-local';EVID=ROOT/'docs/validation/v138/runtime'

def main():
    assert not DEST.exists();DEST.mkdir(parents=True);EVID.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();url='https://codeload.github.com/MmgTools/mmg/tar.gz/refs/tags/v5.8.0'
    op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}));archive=DEST/'mmg-5.8.0.tar.gz'
    with op.open(url,timeout=45) as r,archive.open('wb') as f:
        while True:
            b=r.read(1024**2)
            if not b:break
            f.write(b);assert f.tell()<40*1024**2 and time.monotonic()-start<150
    with tarfile.open(str(archive)) as a:
        for m in a.getmembers():
            if not m.isfile():continue
            parts=Path(m.name).parts;target=(DEST/'source'/Path(*parts[1:])).resolve();target.relative_to((DEST/'source').resolve())
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(a.extractfile(m).read())
    assert sha(DEST/'source/src/mmg3d/libmmg3d.h')==sha(ROOT/'docs/validation/v137/official/mmg-5.8.0-libmmg3d.h')
    files={p.relative_to(DEST).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in DEST.rglob('*') if p.is_file()}
    record=dict(version='5.8.0',url=url,archiveSHA256=sha(archive),publisherIndependentHash=False,files=files,
        seconds=time.monotonic()-start,fileLengthBytes=sum(v['bytes'] for v in files.values()),sourceSHA256=sha(Path(__file__).resolve()),
        status='source-prepared',physicsQualified=False)
    assert record['seconds']<180 and record['fileLengthBytes']<350*1024**2
    (DEST/'source.json').write_text(json.dumps(record,indent=2)+'\n');(EVID/'source.json').write_bytes((DEST/'source.json').read_bytes())
    print(json.dumps({k:v for k,v in record.items() if k!='files'},indent=2))
if __name__=='__main__':main()
