"""独立安装固定CGAL几何审核运行时；不修改系统Python，不上传本地二进制。"""
from pathlib import Path
import hashlib,json,urllib.request,zipfile,subprocess,time,sys
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'tools/runtime/cgal/6.0.1-py311';VERSION='6.0.1.post202410241521'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def get(url):
    with urllib.request.urlopen(url,timeout=45) as r:return r.read()
def extract(file,target):
    with zipfile.ZipFile(str(file)) as z:
        for name in z.namelist():
            p=(target/name).resolve();p.relative_to(target.resolve())
        z.extractall(str(target))
def main():
    assert not DEST.exists();DEST.mkdir(parents=True);start=time.monotonic();sources=[]
    pyurl='https://www.python.org/ftp/python/3.11.9/python-3.11.9-embed-amd64.zip'
    data=json.loads(get('https://pypi.org/pypi/cgal/'+VERSION+'/json'))
    item=next(v for v in data['urls'] if v['filename']=='cgal-'+VERSION+'-cp311-cp311-win_amd64.whl')
    for name,url,expected in [('python.zip',pyurl,None),(item['filename'],item['url'],item['digests']['sha256'])]:
        p=DEST/name;p.write_bytes(get(url));digest=sha(p)
        if expected:assert digest==expected
        sources.append(dict(file=name,url=url,sha256=digest,bytes=p.stat().st_size,matchedPyPIHash=expected is not None))
        extract(p,DEST/'python' if name=='python.zip' else DEST/'packages')
    # 嵌入式运行时只读本地标准库和项目固定包路径，不引入系统site-packages。
    path=DEST/'python/python311._pth';path.write_text('python311.zip\n.\n../packages\n',encoding='ascii')
    command=[str(DEST/'python/python.exe'),'-c','from CGAL import CGAL_Polygon_mesh_processing as p; print(p.self_intersections)']
    proc=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30,creationflags=subprocess.CREATE_NO_WINDOW)
    (DEST/'import.log').write_bytes(proc.stdout)
    files={p.relative_to(DEST).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in DEST.rglob('*') if p.is_file()}
    record=dict(version=VERSION,pythonVersion='3.11.9',exitCode=proc.returncode,seconds=time.monotonic()-start,sources=sources,
        files=files,fileLengthBytes=sum(r['bytes'] for r in files.values()),licenseLocations=[n for n in files if 'license' in n.lower() or 'copying' in n.lower()],
        sourceSHA256=sha(Path(__file__).resolve()))
    (DEST/'runtime.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(proc.stdout.decode('utf8',errors='replace'));print(json.dumps({k:record[k] for k in ['exitCode','seconds','fileLengthBytes','licenseLocations']}))
    assert proc.returncode==0 and record['fileLengthBytes']<=160*1024**2
if __name__=='__main__':main()
