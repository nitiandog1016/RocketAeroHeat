"""复用已转换的原生对象，只编译诊断入口并链接；原引擎安装/核心源码不动。"""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import run,PREFIX,sha
base=ROOT/'outputs/v182_validation';out=base/'build01';start=time.monotonic()
prior=json.loads((out/'reuse_report.json').read_text());assert prior['passed']
used=prior['elapsedSeconds'];records=[]
source=ROOT/'scripts/native/eilmerOperatorDiagnostic.d';compiler=PREFIX/'ldc2-1.41.0-linux-x86_64/bin/ldc2'
src=PREFIX/'source/src';entry_source=out/'executed_diagnostic.d';assert not entry_source.exists();entry_source.write_bytes(source.read_bytes())
for stem,complex_type in [('lmr-run',False),('lmrZ-run',True)]:
    obj=out/(stem+'-native.o');entry=out/(stem+'-entry.o');exe=out/(stem+'-diagnostic')
    flags=['-c','-w','-g','-O2','-release','-boundscheck=off','-dip1008','--d-version=flavour_fast']
    flags += ['--d-version='+x for x in ['enable_fp_exceptions','multi_species_gas','multi_T_gas','MHD','newton_krylov','turbulence']]
    if complex_type:flags+=['--d-version=complex_numbers']
    flags += ['-I'+str(p) for p in [src,src/'extern/D-YAML/source',src/'extern/tinyendian/source',src/'extern/gzip']]
    flags += ['-J='+str(src/'lmr'),'-of='+str(entry)]
    commands=[['objcopy','--redefine-sym=main=rah_original_cmain','--redefine-sym=_Dmain=rah_original_dmain',str(obj)],
              [str(compiler),*flags,str(source)],
              [str(compiler),str(entry),str(obj),str(src/'extern/eqc/libeqc.a'),str(PREFIX/'source/extern/lua-5.4.3/install/lib/liblua.a'),'-L-ldl','-L-lz','-L-Wl,-E','-of='+str(exe)]]
    for i,argv in enumerate(commands):
        remaining=180-used-(time.monotonic()-start);assert remaining>0
        r=run(argv,out,base/(stem+'-build-'+str(i)),min(90,remaining),256,[source,entry_source,src/'lmr'/(stem+'.o'),Path(__file__)])
        records.append(r)
        if not r['passed']:break
    if not records[-1]['passed']:break
report=dict(passed=len(records)==6 and all(r['passed'] for r in records),elapsedSeconds=time.monotonic()-start,
            totalBuildSeconds=used+time.monotonic()-start,processes=records,productionKernelSourceCompiles=0,
            originalEngineObjectsAndInstallPreserved=True,entrySHA256=sha(source))
(out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['passed','elapsedSeconds','totalBuildSeconds','productionKernelSourceCompiles']}))
