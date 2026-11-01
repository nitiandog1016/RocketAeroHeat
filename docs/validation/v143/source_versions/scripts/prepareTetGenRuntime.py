"""下载固定TetGen官方源，并与已取得的Git blob身份逐字核对；不执行网格计算。"""
from pathlib import Path
import hashlib,json,tarfile,time,urllib.request
from prepareCGALRemeshingRuntime import ROOT,sha
DEST=ROOT/'tools/runtime/tetgen/1.5.1';EVID=ROOT/'docs/validation/v142/runtime'

def main():
    start=time.monotonic();assert not DEST.exists();DEST.mkdir(parents=True);EVID.mkdir(parents=True,exist_ok=True)
    metadata=ROOT/'docs/validation/v141/tetgen_research/repository_files.json';rows=json.loads(metadata.read_text())['files']
    url='https://codeload.github.com/TetGen/TetGen/tar.gz/refs/tags/v1.5.1';archive=DEST/'tetgen-1.5.1.tar.gz'
    op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
    with op.open(url,timeout=45) as response,archive.open('wb') as out:
        while True:
            block=response.read(1024**2)
            if not block:break
            out.write(block);assert out.tell()<8*1024**2 and time.monotonic()-start<120
    with tarfile.open(str(archive)) as tar:
        for member in tar.getmembers():
            if not member.isfile():continue
            parts=Path(member.name).parts;target=(DEST/'source'/Path(*parts[1:])).resolve();target.relative_to((DEST/'source').resolve())
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(tar.extractfile(member).read())
    checks={}
    for row in rows:
        p=DEST/'source'/row['name'];data=p.read_bytes();blob=hashlib.sha1(('blob %d\0'%len(data)).encode()+data).hexdigest()
        assert len(data)==row['bytes'] and blob==row['gitBlobSHA1'];checks[row['name']]=dict(sha256=sha(p),gitBlobSHA1=blob,bytes=len(data))
    report=dict(version='1.5.1',url=url,archiveSHA256=sha(archive),publisherArchiveChecksum=False,sourceGitBlobChecks=checks,
        metadataSHA256=sha(metadata),preparationSourceSHA256=sha(Path(__file__).resolve()),seconds=time.monotonic()-start,
        status='prepared-not-compiled-or-tested',physicsQualified=False)
    assert report['seconds']<=120
    (DEST/'source.json').write_text(json.dumps(report,indent=2)+'\n');(EVID/'source.json').write_bytes((DEST/'source.json').read_bytes())
    print(report['seconds'],len(checks),report['archiveSHA256'])

if __name__=='__main__':main()
