"""用真实Windows文件共享锁验证回执原子替换重试和终态备用文件，不启动求解。"""
from pathlib import Path
import ctypes
from ctypes import wintypes
import hashlib
import importlib.util
import json
import threading
import argparse

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='outputs/v112_validation/receipt_lock01');args=parser.parse_args()
    script=ROOT/'scripts/runEmbeddedSU2.py'
    spec=importlib.util.spec_from_file_location('supervisor',str(script));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    out=(ROOT/args.out).resolve();out.relative_to(ROOT/'outputs');out.mkdir(exist_ok=False)
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.CreateFileW.argtypes=[wintypes.LPCWSTR,wintypes.DWORD,wintypes.DWORD,ctypes.c_void_p,wintypes.DWORD,wintypes.DWORD,wintypes.HANDLE]
    kernel.CreateFileW.restype=wintypes.HANDLE;kernel.CloseHandle.argtypes=[wintypes.HANDLE]
    def lock(path):
        h=kernel.CreateFileW(str(path),0x80000000,3,None,3,0,None)
        assert h not in (None,ctypes.c_void_p(-1).value),ctypes.get_last_error()
        return h
    p=out/'transient.json';m.atomic_write_json(p,dict(status='running',step=0));handle=lock(p)
    timer=threading.Timer(.12,lambda:kernel.CloseHandle(handle));timer.start()
    count=m.atomic_write_json(p,dict(status='running',step=1));timer.join()
    assert count>0 and json.loads(p.read_text(encoding='utf-8'))['step']==1
    p=out/'engine_process.json';m.atomic_write_json(p,dict(status='running'));handle=lock(p)
    try:
        failed=False
        try:m.atomic_write_json(p,dict(status='should-not-replace'),maximum_wait=.06)
        except PermissionError:failed=True
        assert failed and json.loads(p.read_text(encoding='utf-8'))['status']=='running'
        fallback=m.save_terminal_receipt(p,dict(status='terminal',exitCode=1,supervisorError='synthetic-write-lock'))
        assert fallback.name=='engine_terminal_fallback.json' and json.loads(fallback.read_text(encoding='utf-8'))['status']=='terminal'
        assert json.loads(p.read_text(encoding='utf-8'))['status']=='running'
    finally:kernel.CloseHandle(handle)
    report=dict(actualWindowsSharingLock=True,transientRetries=count,partialJSONNeverExposed=True,
                persistentLockRaised=True,terminalFallbackPassed=True,solverInvoked=False,
                sourceSHA256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (script,Path(__file__).resolve())})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':main()
