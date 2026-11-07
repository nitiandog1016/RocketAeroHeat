"""用准确编译器布局和原生寄存器读取几何；不修改流动状态。"""
from pathlib import Path
import json
import re
import shutil
import time
from checkEilmerCellTimeScale import ROOT,BASE,PACKAGE,prepare_case
from runEilmerProcess import run,sha

GDB='''set pagination off
set confirm off
set debuginfod enabled off
break _D3lmr11fluidfvcell11FluidFVCell16signal_frequencyMFNibZd
run --max-cpus=1 --max-wall-clock=10
python
import gdb,json,struct,math
off=json.load(open('offsets.json'))
mem=gdb.selected_inferior()
def unpack(addr,fmt): return struct.unpack(fmt,mem.read_memory(addr,struct.calcsize(fmt)))
def value(x): return x if math.isfinite(x) else str(x)
def vector(addr): return list(unpack(addr,'<3d'))
c=int(gdb.parse_and_eval('$rdi'))
n,p=unpack(c+off['lengths'],'<2Q')
assert n==3
fs=unpack(c+off['fs'],'<Q')[0]
nf,pf=unpack(c+off['iface'],'<2Q')
result=dict(observed=True,abi=off,cellAddress=c,
    lengths=[value(x) for x in unpack(p,'<3d')],Lmin=value(unpack(c+off['L_min'],'<d')[0]),
    velocity=vector(fs+off['vel']),soundSpeed=unpack(fs+off['gas']+off['a'],'<d')[0],
    normals=[vector(unpack(pf+8*i,'<Q')[0]+off['normal']) for i in range(nf)],selectedFaces=[0,1,nf-1])
open('time-observation.json','w').write(json.dumps(result,indent=2))
print('RAH-OBSERVATION',json.dumps(result))
end
finish
python
result['signal']=value(float(gdb.parse_and_eval('$xmm0.v2_double[0]')))
open('time-observation.json','w').write(json.dumps(result,indent=2))
print('RAH-SIGNAL',result['signal'])
end
quit
'''

def main():
    """模块一：编译期查询准确ABI，不构建或链接流动内核。"""
    start=time.monotonic();out=BASE/'observation03';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    old=json.loads((ROOT/'outputs/v187_validation/scope03/report.json').read_text())
    fields=dict(lengths='FluidFVCell.lengths',L_min='FluidFVCell.L_min',fs='FluidFVCell.fs',iface='FluidFVCell.iface',
                vel='FlowState.vel',gas='FlowState.gas',a='GasState.a',normal='FVInterface.n')
    source=out/'offsets.d'
    source.write_text('import lmr.fluidfvcell; import lmr.flowstate; import lmr.fvinterface; import gas;\n'+
        '\n'.join('pragma(msg, "RAH-OFFSET '+k+'=", '+v+'.offsetof);' for k,v in fields.items()))
    flags=[x for x in old['compileCommand'][1:] if x.startswith(('-d-version=','--d-version=','-I','-J=')) and x!='-d-version=run_main']
    flags+=['-I/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/extern/'+s for s in ['eqc/source','D-YAML/source','gzip','tinyendian/source']]
    p=run([old['compileCommand'][0],'-c','-o-',*flags,str(source)],out,out/'offsets01',10,32,[source,*old['originalSHA256']],4096,4096)
    assert p['passed']
    offsets={k:int(v) for k,v in re.findall(r'RAH-OFFSET (\w+)=(\d+)',(out/'offsets01/stdout.txt').read_text())}
    assert set(offsets)==set(fields)
    (out/'offsets.json').write_text(json.dumps(offsets,indent=2)+'\n')
    results=[]
    for kind in ['tet','hex']:
        # 模块二：失败控制复用准备；原输入和失败目录始终保持。
        if kind=='tet':
            case=BASE/'tet02';shutil.copytree(BASE/'tet01',case)
        else:
            case=prepare_case(kind)
            for i,cmd in enumerate([['lmr','prep-grids','--job=grid.lua'],['lmr','prep-sim','--job=job.lua']]):
                p=run(cmd,case,out/('hex_prep'+str(i)),min(25,110-(time.monotonic()-start)),32,[Path(__file__),*case.glob('*.*')],4096,4096)
                assert p['passed']
        shutil.copyfile(out/'offsets.json',case/'offsets.json')
        (case/'observe.gdb').write_text(GDB)
        p=run(['gdb','--batch','-x','observe.gdb',str(PACKAGE/'lmr-run')],case,out/(kind+'_native'),
            min(20,110-(time.monotonic()-start)),32,[Path(__file__),PACKAGE/'runtime.json',PACKAGE/'lmr-run',*case.glob('*.*')],4096,4096)
        assert p['passed']
        obs=json.loads((case/'time-observation.json').read_text())
        assert obs['observed'] and 'signal' in obs
        results.append(dict(kind=kind,observation=obs,process=p))
    tet,hexa=[r['observation'] for r in results]
    passed=(tet['lengths']==['nan']*3 and tet['Lmin']>0 and tet['signal']==1e-99 and
            hexa['lengths']==[1.]*3 and hexa['signal']>0)
    report=dict(passed=passed,results=results,elapsedSeconds=time.monotonic()-start,physicalQualification=False,
                nativeUpdates=0,scope='original native object / original first signal; debugger stops inferior')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='results'}))
    return passed

if __name__=='__main__':raise SystemExit(0 if main() else 1)
