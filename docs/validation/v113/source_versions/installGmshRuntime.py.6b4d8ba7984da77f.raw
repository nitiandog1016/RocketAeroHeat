"""安装经SHA固定的Gmsh Python/DLL到项目本地，保留许可证，不改系统环境。"""
from pathlib import Path
import hashlib
import json
import urllib.request
import zipfile

ROOT=Path(__file__).resolve().parents[1]
VERSION='4.15.2'
URL='https://files.pythonhosted.org/packages/78/e9/63c763c0d75544680800651d5ea8c7c3203ae4b5e054706f19a54021a7fb/gmsh-4.15.2-py2.py3-none-win_amd64.whl'
DIGEST='7b36083bb410fa27c5d0e052929d1a9844a5b09169d66017b72b41aabd49d711'
SIZE=42235999

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=ROOT/'tools/runtime/gmsh'/VERSION;target.mkdir(parents=True,exist_ok=True)
    record=target/'runtime.json'
    if record.exists():
        m=json.loads(record.read_text(encoding='utf-8'))
        assert m['archiveSHA256']==DIGEST
        for f in m['files']:assert sha(target/f['path'])==f['sha256']
        print('Existing Gmsh runtime verified');return
    package=target/'official.whl'
    if not package.exists():
        partial=target/'download.part'
        with urllib.request.urlopen(URL,timeout=45) as src,partial.open('wb') as dst:
            count=0
            while True:
                b=src.read(1024**2)
                if not b:break
                count+=len(b);assert count<=SIZE;dst.write(b)
        assert partial.stat().st_size==SIZE and sha(partial)==DIGEST
        partial.replace(package)
    assert sha(package)==DIGEST
    files=[]
    with zipfile.ZipFile(str(package)) as z:
        selected=[i for i in z.infolist() if not i.is_dir() and
                  (Path(i.filename).name=='gmsh.py' or i.filename.lower().endswith('.dll') or
                   Path(i.filename).name.upper().startswith(('LICENSE','COPYING')))]
        assert sum(i.file_size for i in selected)<200*1024**2
        for i in selected:
            p=target/Path(i.filename).name;b=z.read(i)
            if p.exists():assert p.read_bytes()==b
            else:p.write_bytes(b)
            files.append(dict(path=p.name,sha256=sha(p),bytes=len(b)))
    assert (target/'gmsh.py').is_file() and any(f['path'].endswith('.dll') for f in files)
    m=dict(version=VERSION,archiveURL=URL,archiveSHA256=DIGEST,archiveBytes=SIZE,
           files=files,license='GPL-2.0-or-later; upstream exceptions apply',
           sourceURL='https://gmsh.info/',systemPATHModified=False,physicsQualified=False)
    record.write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(m))

if __name__=='__main__':main()
