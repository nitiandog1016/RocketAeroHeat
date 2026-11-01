"""为局部体网格研究准备固定C++头文件；不修改系统或旧Python自交运行时。"""
from pathlib import Path
import hashlib,json,os,subprocess,tarfile,time,urllib.request
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'tools/runtime/cgal/6.0.1-native'
EVID=ROOT/'docs/validation/v137/runtime'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    assert not DEST.exists();DEST.mkdir(parents=True);EVID.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();record=dict(status='preparing',sources=[],physicsQualified=False)
    op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
    try:
        for name,url in [('CGAL-6.0.1.tar.xz','https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.tar.xz'),
                         ('boost_1_87_0.tar.bz2','https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.bz2')]:
            p=DEST/name
            with op.open(url,timeout=45) as r,p.open('wb') as f:
                while True:
                    b=r.read(1024*1024)
                    if not b:break
                    f.write(b)
                    assert p.stat().st_size<160*1024**2 and time.monotonic()-start<450
            item=dict(file=name,url=url,sha256=sha(p),bytes=p.stat().st_size,independentPublisherHash=False)
            if name.startswith('boost'):
                meta=op.open(url+'.json',timeout=30).read();(EVID/'boost_archive.json').write_bytes(meta)
                expected=json.loads(meta)['sha256'];assert item['sha256']==expected
                item['independentPublisherHash']=True
            record['sources'].append(item)
            with tarfile.open(str(p)) as archive:
                for m in archive.getmembers():
                    parts=Path(m.name).parts
                    if len(parts)<2:continue
                    rel=Path(*parts[1:])
                    keep=(parts[1] in ['include','boost','cmake']) or ('LICENSE' in m.name or 'COPYING' in m.name)
                    if not keep or not m.isfile():continue
                    target=(DEST/('cgal' if name.startswith('CGAL') else 'boost')/rel).resolve()
                    target.relative_to(DEST.resolve());target.parent.mkdir(parents=True,exist_ok=True)
                    target.write_bytes(archive.extractfile(m).read())
            print('prepared '+name,flush=True)
        compiler=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe'
        record.update(status='headers-prepared',compiler=str(compiler),compilerSHA256=sha(compiler),
            compilerVersion=subprocess.check_output([str(compiler),'--version'],creationflags=subprocess.CREATE_NO_WINDOW).decode(),
            sourceSHA256=sha(Path(__file__)),seconds=time.monotonic()-start)
        files={p.relative_to(DEST).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in DEST.rglob('*') if p.is_file()}
        record['files']=files;record['fileLengthBytes']=sum(v['bytes'] for v in files.values())
        assert record['fileLengthBytes']<650*1024**2
    except Exception as e:
        record.update(status='failed',error=repr(e),seconds=time.monotonic()-start)
        raise
    finally:
        (DEST/'headers.json').write_text(json.dumps(record,indent=2)+'\n')
        (EVID/'headers.json').write_bytes((DEST/'headers.json').read_bytes())
if __name__=='__main__':main()
