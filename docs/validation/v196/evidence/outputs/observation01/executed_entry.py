"""原生单元几何时间尺度核验；只观察第一次信号，不推进整箭。"""
from pathlib import Path
import json
import shutil
import sys
import time
import numpy as np
from runEilmerProcess import run, sha

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v196_validation'
SOURCE=ROOT/'outputs/v194_validation/rocket01'
PACKAGE=ROOT/'tools/runtime/eilmer-explicit-workspace-v1'

GDB='''set pagination off
set confirm off
set print elements 8
set print max-depth 4
set debuginfod enabled off
break _D3lmr11fluidfvcell11FluidFVCell16signal_frequencyMFNibZd
run --max-cpus=1 --max-wall-clock=10
info args
python
import gdb,json,math
result={}
try:
    c=gdb.selected_frame().read_var('this').dereference()
    def scalar(x):
        y=float(x)
        return y if math.isfinite(y) else str(y)
    def vec(x): return [scalar(x[k]) for k in ('x','y','z')]
    def arr(x):
        names=[f.name for f in x.type.fields()]
        result.setdefault('sliceFields',names)
        return x['ptr'],int(x['length'])
    p,n=arr(c['lengths'])
    result['lengths']=[scalar(p[i]) for i in range(n)]
    result['Lmin']=scalar(c['L_min'])
    fs=c['fs'].dereference()
    result['velocity']=vec(fs['vel'])
    result['soundSpeed']=scalar(fs['gas']['a'])
    p,n=arr(c['iface'])
    result['normals']=[vec(p[i].dereference()['n']) for i in range(n)]
    result['selectedFaces']=[0,1,n-1]
    result['observed']=True
except Exception as e:
    result['error']=repr(e)
    print(gdb.execute('print *this',to_string=True))
open('time-observation.json','w').write(json.dumps(result,indent=2))
print('RAH-OBSERVATION',json.dumps(result))
end
finish
python
try:
    result['signal']=scalar(gdb.parse_and_eval('$xmm0.v2_double[0]'))
    open('time-observation.json','w').write(json.dumps(result,indent=2))
    print('RAH-SIGNAL',result['signal'])
except Exception as e: print('RAH-SIGNAL-ERROR',repr(e))
end
quit
'''

def prepare_case(kind):
    """模块一：一颗真实原四面体及独立六面体，全部环境边界保持均匀状态。"""
    case=BASE/(kind+'01');case.mkdir(parents=True,exist_ok=False)
    if kind=='tet':
        z=np.load(ROOT/'outputs/v194_validation/netgen_recovered02/mesh.npz')
        ids=z['tetrahedra'][0];v=z['vertices'][ids]
        typ=10;faces=[[0,2,1],[0,1,3],[1,2,3],[2,0,3]];ft=5
        origin=dict(meshSHA256=sha(ROOT/'outputs/v194_validation/netgen_recovered02/mesh.npz'),cellIndex=0,vertexIDs=ids.tolist())
    else:
        v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
        typ=12;faces=[[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7],[0,3,2,1],[4,5,6,7]];ft=9
        origin=dict(geometry='unit cube; independent of rocket')
    text=['NDIME= 3','NELEM= 1',str(typ)+' '+' '.join(map(str,range(len(v))))+' 0','NPOIN= '+str(len(v))]
    text+=[' '.join(format(x,'.17g') for x in row)+' '+str(i) for i,row in enumerate(v)]
    text+=['NMARK= 1','MARKER_TAG= farfield','MARKER_ELEMS= '+str(len(faces))]
    text+=[str(ft)+' '+' '.join(map(str,f)) for f in faces]
    (case/'mesh.su2').write_text('\n'.join(text)+'\n')
    for name in ['rah-air.gas','grid.lua']:
        shutil.copyfile(SOURCE/name,case/name)
    job=(SOURCE/'job.lua').read_text()
    job=job.replace('wall=WallBC_WithSlip:new{group="wall"},\n farfield=', 'farfield=')
    job=job.replace('config.max_step=0','config.max_step=1')
    (case/'job.lua').write_text(job)
    (case/'origin.json').write_text(json.dumps(origin,indent=2)+'\n')
    (case/'observe.gdb').write_text(GDB)
    return case

def main():
    """模块二：固定总预算、准确源码与输入锁定；断点退出只记观察完成。"""
    started=time.monotonic();BASE.mkdir(exist_ok=True)
    out=BASE/'observation01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    package=json.loads((PACKAGE/'runtime.json').read_text())
    assert all(sha(p)==h for p,h in package['filesSHA256'].items())
    results=[]
    for kind in ['tet','hex']:
        case=prepare_case(kind)
        for index,cmd in enumerate([['lmr','prep-grids','--job=grid.lua'],['lmr','prep-sim','--job=job.lua']]):
            left=min(25,120-(time.monotonic()-started));assert left>0
            p=run(cmd,case,out/(kind+'_prep'+str(index)),left,32,[Path(__file__),*case.glob('*.*')],4096,4096)
            assert p['passed'],p.get('stopReason')
        p=run(['gdb','--batch','-x','observe.gdb',str(PACKAGE/'lmr-run')],case,out/(kind+'_native'),
              min(20,120-(time.monotonic()-started)),32,
              [Path(__file__),PACKAGE/'runtime.json',PACKAGE/'lmr-run',*case.glob('*.*')],4096,4096)
        assert p['passed'],'调试过程未按预算结束'
        obs=json.loads((case/'time-observation.json').read_text())
        results.append(dict(kind=kind,observation=obs,process=p))
    report=dict(passed=all(x['observation'].get('observed',False) and 'signal' in x['observation'] for x in results),
        results=results,elapsedSeconds=time.monotonic()-started,physicalQualification=False,nativeUpdates=0,
        scope='First native signal; debugger quits stopped inferior, not a successful flow solve')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='results'}))
    return report['passed']

if __name__=='__main__':raise SystemExit(0 if main() else 1)
