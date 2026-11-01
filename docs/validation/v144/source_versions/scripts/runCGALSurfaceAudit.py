"""MATLAB/Python主后端可调用的有界CGAL自交检查；实际进程终态与原始输入锁定。"""
from pathlib import Path
import json,os,subprocess,time
import numpy as np
from scipy.io import savemat
from prepareSU2CheckpointContinuation import ROOT,read,sha
RUNTIME=ROOT/'tools/runtime/cgal/6.0.1-py311'

def audit_surface(P,F,folder,source_metadata):
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs');assert not folder.exists();folder.mkdir(parents=True)
    P=np.asarray(P);F=np.asarray(F,dtype=int);assert P.ndim==2 and P.shape[1]==3 and np.isfinite(P).all()
    assert F.shape[1]==3 and np.min(F)>=0 and np.max(F)<len(P)
    used,inverse=np.unique(F,return_inverse=True);P=P[used];F=inverse.reshape(F.shape)
    assert np.all(np.linalg.norm(np.cross(P[F[:,1]]-P[F[:,0]],P[F[:,2]]-P[F[:,0]]),axis=1)>0)
    savemat(str(folder/'audit_input.mat'),dict(V=P,F=F+1,sourceNodeIDsZeroBased=used),do_compression=True)
    off=folder/'surface.off'
    with off.open('w',encoding='ascii',newline='\n') as stream:
        stream.write('OFF\n%d %d 0\n'%(len(P),len(F)));np.savetxt(stream,P,fmt='%.17g');np.savetxt(stream,np.c_[np.full(len(F),3),F],fmt='%d')
    runtime=read(RUNTIME/'runtime.json');assert runtime['exitCode']==0
    for name,row in runtime['files'].items():assert sha(RUNTIME/name)==row['sha256']
    files=[off,folder/'audit_input.mat',Path(__file__).resolve(),ROOT/'scripts/cgalSurfaceIntersections.py',RUNTIME/'runtime.json']
    before={str(p):sha(p) for p in files};receipt=dict(status='running',sourceSHA256=before,sourceMetadata=source_metadata,maximumSeconds=60,supervisorPID=os.getpid())
    target=folder/'process.json'
    def save():target.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    save();start=time.monotonic();timed=False
    with (folder/'stdout.log').open('wb') as log:
        process=subprocess.Popen([str(RUNTIME/'python/python.exe'),str(ROOT/'scripts/cgalSurfaceIntersections.py'),'--input',str(off),'--output',str(folder/'intersections.json')],
            stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
        receipt['processID']=process.pid;save()
        try:code=process.wait(timeout=60)
        except subprocess.TimeoutExpired:timed=True;process.kill();code=process.wait()
    stable={name:sha(Path(name))==h for name,h in before.items()}
    runtimeStable=all(sha(RUNTIME/name)==row['sha256'] for name,row in runtime['files'].items())
    receipt.update(status='terminal',exitCode=code,timedOut=timed,seconds=time.monotonic()-start,sourceStableAfter=stable,
        runtimeStableAfter=runtimeStable,outputSHA256={p.name:sha(p) for p in folder.iterdir() if p.is_file() and p!=target})
    save();assert code==0 and not timed and all(stable.values()) and runtimeStable
    result=read(folder/'intersections.json');assert result['inputSHA256']==sha(off)
    return result
