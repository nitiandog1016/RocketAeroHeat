"""定向绑定复数NK/SA分配源并查询编译期布局，不构建或执行流动引擎。"""
from pathlib import Path
import json,re,shutil,time
from runEilmerProcess import run,sha,PREFIX,COMMIT

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v211_validation'

def save(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def main():
    """模块一：只提取准确编译来源；保留实际执行入口，供离线账复核。"""
    start=time.monotonic();out=BASE/'source01';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    prior=ROOT/'outputs/v197_validation/build01/report.json'
    r=json.loads(prior.read_text());assert r['passed'] and r['inputsStable']
    src=PREFIX/'source/src'
    names=['lmr/'+s+'.d' for s in ['fluidblock','ufluidblock','coredata','fluidfvcell','fvinterface',
        'fvvertex','flowstate','flowgradients','lsqinterp','onedinterp','newtonkrylovsolver',
        'conservedquantities','globalconfig','turbulence/turb_model']]
    # 所有布局相关D源均从原成功构建的清单定位，不猜测湍流子目录。
    names=[n for n in names if (src/n).is_file()]
    names+=['nm/bbla.d','nm/smla.d','nm/number.d','geom/grid/usgrid.d','gas/gas_state.d']
    selected=[]
    for name in names:
        p=src/name
        if name=='lmr/fluidfvcell.d':p=ROOT/'tools/runtime/eilmer-nk-geometric-v1/overlay/src/lmr/fluidfvcell.d'
        assert sha(p)==r['inputSHA256'][str(p)], str(p)
        q=out/'source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        selected.append(dict(path=str(p),relativePath=name,sha256=sha(p)))
    # 已有构建清单中湍流配置/构造模块，保存供六变量及分配路径核对。
    for name,h in r['inputSHA256'].items():
        p=Path(name)
        if '/lmr/' in name and p.suffix=='.d' and ('turb' in p.name or p.name=='blockconfig.d'):
            assert sha(p)==h
            q=out/'source'/p.relative_to(src);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
            selected.append(dict(path=name,relativePath=p.relative_to(src).as_posix(),sha256=h))
    source_report=dict(passed=True,sourceCommit=COMMIT,files=selected,seconds=time.monotonic()-start,
        engineBuilds=0,CFDRuns=0)
    assert source_report['seconds']<=60
    save(out/'report.json',source_report)
    # 模块二：只编译类型查询，无对象文件、无链接，也不实例化网格或求解器。
    query=BASE/'layout01';query.mkdir(exist_ok=False)
    types={
        'FlowGradients':'FlowGradients.sizeof','LSQInterpGradients':'LSQInterpGradients.sizeof',
        'WLSQGradWorkspace':'WLSQGradWorkspace.sizeof','LSQInterpWorkspace':'LSQInterpWorkspace.sizeof',
        'FlowState':'FlowState.sizeof','FluidFVCell':'__traits(classInstanceSize, FluidFVCell)',
        'FVInterface':'__traits(classInstanceSize, FVInterface)','FVVertex':'__traits(classInstanceSize, FVVertex)',
        'USGFace':'__traits(classInstanceSize, USGFace)','USGCell':'__traits(classInstanceSize, USGCell)',
        'Matrix':'__traits(classInstanceSize, Matrix!double)',
        'PhaseAir':'__traits(classInstanceSize, RAHPhaseAwareAir)',
        'Vector3':'Vector3.sizeof',
        'Slice':'(number[]).sizeof','Scalar':'number.sizeof'}
    probe=query/'layout.d'
    probe.write_text('module rah_nksa_layout;\n'
        'import lmr.flowgradients; import lmr.lsqinterp; import lmr.flowstate;\n'
        'import lmr.fluidfvcell; import lmr.fvinterface; import lmr.fvvertex; import geom.grid.usgrid;\n'
        'import nm.bbla; import nm.number; import geom; import gas.rah_phase_aware_air;\n'+
        '\n'.join('pragma(msg, "RAH-LAYOUT %s=", %s);'%(k,v) for k,v in types.items())+'\n')
    flags=[s for s in r['argv'][1:] if s.startswith(('-d-version=','--d-version=','-I','-J')) and s!='-d-version=run_main']
    overlay=ROOT/'outputs/v210_validation/prepared02/overlay'
    flags=['-I'+str(overlay),'-I'+str(ROOT/'tools/runtime/eilmer-nk-geometric-v1/overlay/src'),*flags]
    flags+=['-I'+str(src/'extern'/p) for p in ['eqc/source','D-YAML/source','gzip','tinyendian/source']]
    flags+=['-J'+str(ROOT/'data/thermo')]
    command=[r['argv'][0],'-c','-o-','-boundscheck=on',*flags,str(probe)]
    assert all(f in command for f in ['-d-version=complex_numbers','-d-version=newton_krylov','-d-version=multi_species_gas'])
    locks=[Path(p) for p in r['inputSHA256'] if Path(p).is_file()]
    locks+=[p for p in overlay.rglob('*.d') if p.is_file()]
    locks+=[probe,Path(__file__),ROOT/'docs/validation/v211/PLAN.md',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    save(query/'admission.json',dict(seconds=30,parentSeconds=45,memoryMiB=2048,outputMiB=16,automaticRetry=0,command=command))
    receipt=run(command,query,BASE/'layout_process01',30,16,locks+[query/'admission.json'],2048,2048)
    sizes={k:int(v) for k,v in re.findall(r'RAH-LAYOUT (\w+)=(\d+)',(BASE/'layout_process01/stdout.txt').read_text())}
    report=dict(passed=receipt['passed'] and set(sizes)==set(types),sizes=sizes,process=receipt,
        seconds=time.monotonic()-start,engineBuilds=0,CFDRuns=0,physicalQualification=False)
    save(query/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='process'}))
    return report['passed']

if __name__=='__main__':
    raise ValueError('v211定向源和原生布局查询已完成，CLI关闭；不重复查询或构建。')
