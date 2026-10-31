"""固定可压缩CoolProp构建依赖；仅项目内安装，保留来源、许可与原字节。"""
from pathlib import Path
import hashlib,json,urllib.request,zipfile,concurrent.futures,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'tools/runtime/su2-coolprop-build'
OLD=ROOT/'tools/runtime/su2-thermal-build';COOL='98b3523d5daa98454618d381d2ae53f7471d216b'
DEPS=[('Eigen','libeigen/eigen','68f4e58cfacc686583d16cff90361f0b43bc2c1b'),
 ('IF97','CoolProp/IF97','7aaced024a702f0985474bf293cdaae9c8d06521'),
 ('REFPROP-headers','CoolProp/REFPROP-headers','b4faab1b73911c32c4b69c526c7e92f74edb67de'),
 ('fmtlib','fmtlib/fmt','e424e3f2e607da02742f73db84873b8084fc714c'),
 ('msgpack-c','msgpack/msgpack-c','919908742b4fdbc575e77fe1a8657e70c9573c44'),
 ('rapidjson','Tencent/rapidjson','24b5e7a8b27f42fa16b96fc70aade9106cf7102f')]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def download(name,url,limit,expected=None):
    file=BASE/(name+'.zip');receipt=BASE/(name+'_download.json')
    if receipt.exists():
        r=json.loads(receipt.read_text());assert sha(file)==r['sha256'];return r
    assert not file.exists() and not file.with_suffix('.partial').exists()
    start=time.monotonic();partial=file.with_suffix('.partial')
    opener=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
    with opener.open(url,timeout=45) as src,partial.open('xb') as dst:
        count=0
        while True:
            b=src.read(1024**2)
            if not b:break
            count+=len(b);assert count<=limit and time.monotonic()-start<=180;dst.write(b)
    digest=sha(partial)
    if expected:assert digest==expected
    partial.replace(file)
    r=dict(url=url,sha256=digest,bytes=count,upstreamDigestChecked=bool(expected),elapsedSeconds=time.monotonic()-start)
    receipt.write_text(json.dumps(r,indent=2)+'\n');return r
def extract(name,target):
    marker=target/'.rah_extracted.json';archive=BASE/(name+'.zip')
    if marker.exists():assert json.loads(marker.read_text())['archiveSHA256']==sha(archive);return
    with zipfile.ZipFile(str(archive)) as z:
        items=z.infolist();top=items[0].filename.split('/')[0];assert sum(q.file_size for q in items)<=600*1024**2
        for q in items:
            if q.is_dir():continue
            rel=Path(q.filename).relative_to(top);p=(target/rel).resolve();p.relative_to(target.resolve())
            assert (q.external_attr>>16)&0o170000!=0o120000,'须单独审核符号链接'
            p.parent.mkdir(parents=True,exist_ok=True);data=z.read(q)
            if p.exists():assert p.read_bytes()==data
            else:p.write_bytes(data)
    marker.write_text(json.dumps(dict(archiveSHA256=sha(archive)))+'\n')
def main():
    BASE.mkdir(parents=True,exist_ok=True)
    opener=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
    sums=opener.open('https://github.com/Kitware/CMake/releases/download/v3.31.8/cmake-3.31.8-SHA-256.txt',timeout=30).read()
    (BASE/'cmake-3.31.8-SHA-256.txt').write_bytes(sums)
    expected=next(line.split()[0] for line in sums.decode().splitlines() if line.endswith('cmake-3.31.8-windows-x86_64.zip'))
    jobs=[('cmake','https://github.com/Kitware/CMake/releases/download/v3.31.8/cmake-3.31.8-windows-x86_64.zip',60*1024**2,expected),
          ('CoolProp','https://github.com/CoolProp/CoolProp/archive/'+COOL+'.zip',100*1024**2,None)]
    for name,repo,commit in DEPS:
        url='https://github.com/'+repo+'/archive/'+commit+'.zip'
        if name=='Eigen':url='https://gitlab.com/libeigen/eigen/-/archive/'+commit+'/eigen-'+commit+'.zip'
        jobs.append((name,url,30*1024**2,None))
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(download,*j) for j in jobs];results=[];errors=[]
        for job,f in zip(jobs,futures):
            try:r=f.result();results.append(r);print(json.dumps(dict(name=job[0],**r)),flush=True)
            except Exception as e:errors.append(dict(name=job[0],error=repr(e)))
    if errors:
        (BASE/'download_failure.json').write_text(json.dumps(errors,indent=2)+'\n');raise RuntimeError(errors)
    assert sum(r['bytes'] for r in results)<=250*1024**2
    extract('cmake',BASE/'cmake');extract('CoolProp',OLD/'source/subprojects/CoolProp')
    for name,_,_ in DEPS:extract(name,OLD/'source/subprojects/CoolProp/externals'/name)
    report=dict(su2Commit='7048fd35f56d4bf52e586429d6746a36a10255f8',coolPropCommit=COOL,cmakeVersion='3.31.8',
        dependencies=[dict(name=n,repository=r,commit=c) for n,r,c in DEPS],downloadBytes=sum(r['bytes'] for r in results),
        projectLocalOnly=True,physicalSolutionQualified=False)
    (BASE/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');print('Fixed build inputs prepared; no engine compiled.',flush=True)
if __name__=='__main__':main()
