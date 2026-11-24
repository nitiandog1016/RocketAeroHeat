"""安装经SHA固定的Gmsh Python/DLL到项目本地，保留许可证，不改系统环境。"""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import http.client
import json
import urllib.request
import zipfile

ROOT=Path(__file__).resolve().parents[1]
VERSION='4.15.2'
URL='https://files.pythonhosted.org/packages/78/e9/63c763c0d75544680800651d5ea8c7c3203ae4b5e054706f19a54021a7fb/gmsh-4.15.2-py2.py3-none-win_amd64.whl'
DIGEST='7b36083bb410fa27c5d0e052929d1a9844a5b09169d66017b72b41aabd49d711'
SIZE=42235999

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def download_pinned(url, package, partial, size, digest):
    """按服务器确认的字节范围续传，只有完整官方摘要相符才发布安装包。"""
    chunk = 2 * 1024**2
    with partial.open('ab+') as dst:
        offset = dst.tell()
        if offset > size:
            raise ValueError(f'暂存文件超过登记大小：{partial}')
        while offset < size:
            end = min(offset + chunk, size) - 1
            request = urllib.request.Request(url, headers={'Range': f'bytes={offset}-{end}', 'Accept-Encoding': 'identity'})
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(request, timeout=45) as src:
                        full = src.status == 200
                        if not full and (src.status != 206 or src.headers.get('Content-Range') != f'bytes {offset}-{end}/{size}'):
                            raise ValueError('服务器响应的字节范围与 文件下载请求不符。')
                        expected = size if full else end - offset + 1
                        data = src.read(expected + 1)
                    if len(data) != expected:
                        raise OSError(f'文件下载片段不完整，期望 {expected} 字节，收到 {len(data)}。')
                    break
                except (OSError, http.client.HTTPException):
                    if attempt == 2:
                        raise
            if full:
                if hashlib.sha256(data).hexdigest() != digest:
                    raise ValueError('服务器返回的完整文件摘要错误。')
                dst.truncate(0)
                offset = 0
            dst.write(data)
            dst.flush()
            offset += len(data)
    actual = sha(partial)
    if actual != digest:
        # 完整但损坏的文件留存。下次从新的暂存文件重试。
        damaged = partial.with_name('download-'+actual+'.invalid')
        partial.replace(damaged)
        raise ValueError(f'文件摘要不符，已保留 {damaged}。')
    partial.replace(package)

def main():
    target=ROOT/'tools/runtime/gmsh'/VERSION;target.mkdir(parents=True,exist_ok=True)
    if (target/'runtime.json').is_file():
        install(target)
        return
    with file_lock(target/'install.lock'):
        install(target)


@contextmanager
def file_lock(path):
    """进程退出时由系统释放文件锁，保护同一暂存包和最终发布文件。"""
    import msvcrt
    with path.open('a+b') as lock:
        lock.seek(0)
        try:
            msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        except OSError as exc:
            raise RuntimeError(f'另一个进程正在准备 {path.parent}。') from exc
        try:
            yield
        finally:
            lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)

def install(target):
    record=target/'runtime.json'
    if record.exists():
        m=json.loads(record.read_text(encoding='utf-8'))
        assert m['archiveSHA256']==DIGEST
        for f in m['files']:assert sha(target/f['path'])==f['sha256']
        print('Existing Gmsh runtime verified');return
    package=target/'official.whl'
    if not package.exists():
        download_pinned(URL, package, target/'download.part', SIZE, DIGEST)
    assert sha(package)==DIGEST
    files=[]
    with zipfile.ZipFile(str(package)) as z:
        selected=[i for i in z.infolist() if not i.is_dir() and
                  (Path(i.filename).name=='gmsh.py' or i.filename.lower().endswith('.dll') or
                   Path(i.filename).name.upper().startswith(('LICENSE','COPYING')))]
        assert sum(i.file_size for i in selected)<200*1024**2
        for i in selected:
            p=target/Path(i.filename).name;b=z.read(i)
            if not p.exists() or p.read_bytes()!=b:
                pending=p.with_name(p.name+'.pending');pending.write_bytes(b);pending.replace(p)
            files.append(dict(path=p.name,sha256=sha(p),bytes=len(b)))
    assert (target/'gmsh.py').is_file() and any(f['path'].endswith('.dll') for f in files)
    m=dict(version=VERSION,archiveURL=URL,archiveSHA256=DIGEST,archiveBytes=SIZE,
           files=files,license='GPL-2.0-or-later; upstream exceptions apply',
           sourceURL='https://gmsh.info/',systemPATHModified=False,physicsQualified=False)
    pending=record.with_suffix('.pending')
    pending.write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8');pending.replace(record)
    print(json.dumps(m))

if __name__=='__main__':main()
