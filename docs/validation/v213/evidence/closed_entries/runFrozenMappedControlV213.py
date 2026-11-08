"""分阶段监督原生接口控制；每个阶段使用独立小cwd，固定累计预算。"""
from pathlib import Path
import json,shutil,sys,time
from runEilmerProcess import run,sha,PREFIX,COMMIT
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v213_validation';PACKAGE=ROOT/'tools/runtime/eilmer-frozen-mapped-probe-v1'
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')

def geometry():
    case=BASE/'control01';case.mkdir(exist_ok=False)
    locks=[ROOT/'scripts/prepareFrozenMappedControlV213.py',ROOT/'docs/validation/v213/PLAN.md',ROOT/'outputs/v208_validation/layer_mesh02/exact_mesh.npz']
    a=run([sys.executable,str(locks[0])],case,BASE/'geometry01',60,32,locks,1024,1024);assert a['passed']
    a=run(['lmr','prep-grids','--job=grid.lua'],case,BASE/'prep_grid01',30,32,
        locks+[case/'mesh.su2',case/'grid.lua'],2048,2048);assert a['passed']

def build():
    """复用已成功的完整清单；只换main/具名气体登记/新边界，原内核不改。"""
    start=time.monotonic();PACKAGE.mkdir(exist_ok=True)
    previous=json.loads((ROOT/'outputs/v197_validation/build01/report.json').read_text())
    assert previous['passed'] and previous['sourceCommit']==COMMIT
    for name,h in previous['inputSHA256'].items():assert sha(name)==h
    shutil.copyfile(ROOT/'tools/runtime/eilmer-nk-geometric-v1/buildinfo.json',PACKAGE/'buildinfo.json')
    native=ROOT/'native/eilmerFrozenMapped';gas=ROOT/'outputs/v210_validation/prepared02/overlay'
    command=[]
    for arg in previous['argv']:
        if arg=='-d-version=run_main':continue
        if arg.startswith('-of='):arg='-of='+str(PACKAGE/'probe')
        if arg.endswith('/src/gas/init_gas_model.d'):arg=str(gas/'gas/init_gas_model.d')
        if arg.startswith('-J=') and 'eilmer-nk-geometric-v1' in arg:arg='-J='+str(PACKAGE)
        command.append(arg)
    command[1:1]=['-I'+str(native),'-I'+str(gas),'-J'+str(ROOT/'data/thermo')]
    command += [str(native/'probe.d'),str(native/'lmr/bc/ghost_cell_effect/rah_frozen_mapped.d'),str(gas/'gas/rah_phase_aware_air.d')]
    locks=list(previous['inputSHA256'])+[p for p in native.rglob('*') if p.is_file()]+[p for p in gas.rglob('*') if p.is_file()]
    locks += [ROOT/'docs/validation/v213/PLAN.md',Path(__file__).resolve(),PACKAGE/'buildinfo.json',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    old=sorted(BASE.glob('build*/report.json'));spent=sum(json.loads(p.read_text())['elapsedSeconds'] for p in old)
    assert spent<240;record=BASE/('build%02d'%(len(old)+1))
    process=run(command,PACKAGE,record,240-spent,512,locks,6144,6144)
    receipt=dict(passed=process['passed'],process=process,entrySeconds=time.monotonic()-start,
        cumulativeBuildSeconds=spent+process['elapsedSeconds'],sourceCommit=COMMIT,CFDAdvanceRuns=0,
        profile='native-frozen-mapped-operator-probe-v1',productionFactoryRegistered=False)
    if process['passed']:receipt['binarySHA256']=sha(PACKAGE/'probe')
    save(BASE/'build_result.json',receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k!='process'}));return process['passed']

def prepare():
    r=json.loads((BASE/'build_result.json').read_text());assert r['passed'] and sha(PACKAGE/'probe')==r['binarySHA256']
    case=BASE/'control01';binary=PACKAGE/'bin/probe';assert sha(binary)==r['binarySHA256']
    locks=[p for p in case.rglob('*') if p.is_file()]+[binary,PACKAGE/'lib/prepsim.lua',
        ROOT/'docs/validation/v213/PLAN.md',ROOT/'docs/validation/v213/PREPARATION_RECOVERY.md']
    spent=json.loads((BASE/'prep_sim01/report.json').read_text())['elapsedSeconds']
    r=run([str(binary),'prepare'],case,BASE/'prep_sim02',30-spent,32,locks,2048,2048)
    print(json.dumps(dict(passed=r['passed'],seconds=r['elapsedSeconds'])));return r['passed']

def query():
    case=BASE/'control01';r=json.loads((BASE/'build_result.json').read_text());assert r['passed'] and sha(PACKAGE/'probe')==r['binarySHA256']
    assert json.loads((BASE/'prep_sim02/report.json').read_text())['passed']
    locks=[p for p in case.rglob('*') if p.is_file()]+[PACKAGE/'bin/probe',ROOT/'docs/validation/v213/PLAN.md']
    r=run([str(PACKAGE/'bin/probe')],case,BASE/'query01',60,32,locks,2048,2048)
    print(json.dumps(dict(passed=r['passed'],seconds=r['elapsedSeconds'],exit=r['nativeExitCode'])));return r['passed']

if __name__=='__main__':
    stage=sys.argv[1]
    if stage=='geometry':geometry()
    elif stage=='build':raise SystemExit(0 if build() else 1)
    elif stage=='prepare':raise SystemExit(0 if prepare() else 1)
    elif stage=='query':raise SystemExit(0 if query() else 1)
    else:raise ValueError('unknown stage')
