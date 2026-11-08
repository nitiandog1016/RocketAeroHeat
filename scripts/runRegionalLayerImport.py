"""监督唯一原区域层转换与checkMesh；不启动任何流动计算。"""
import json,os,sys,time,subprocess,shutil
from pathlib import Path
from prepareOpenFoamControl import ROOT,sha,save
BASE=ROOT/'outputs/v208_validation'

def export():
    import psutil
    out=BASE/'layer_export01';out.mkdir(exist_ok=False)
    paths=[ROOT/'scripts/exportRegionalPrismLayers.py',Path(__file__).resolve(),ROOT/'docs/validation/v208/LAYER_IMPORT_PLAN.md',
        ROOT/'outputs/v144_validation/actual/hybrid_volume.mat',BASE/'layers02/selection.npz']
    before={str(p.relative_to(ROOT)):sha(p) for p in paths}
    save(out/'admission.json',dict(seconds=120,memoryMiB=2048,outputMiB=256,inputSHA256=before,automaticRetry=0))
    shutil.copyfile(__file__,out/'executed_entry.py');start=time.monotonic();peak=0;stop=None
    with (out/'stdout.txt').open('wb') as stream:
        env=dict(os.environ,PYTHONIOENCODING='utf-8',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        child=subprocess.Popen([sys.executable,str(paths[0])],stdout=stream,stderr=subprocess.STDOUT,env=env)
        proc=psutil.Process(child.pid)
        while child.poll() is None:
            try:peak=max(peak,proc.memory_info().rss)
            except psutil.NoSuchProcess:pass
            elapsed=time.monotonic()-start
            size=sum(p.stat().st_size for p in (BASE/'layer_mesh01').rglob('*') if p.is_file()) if (BASE/'layer_mesh01').exists() else 0
            if elapsed>120 or peak>2048*1024**2 or size>256*1024**2:
                stop='resource-limit';child.kill();break
            time.sleep(.1)
        code=child.wait()
    stable=all(sha(ROOT/p)==h for p,h in before.items())
    result=dict(passed=code==0 and stop is None and stable,exitCode=code,stopReason=stop,elapsedSeconds=time.monotonic()-start,
        peakObservedRSSBytes=peak,inputStable=stable,inputSHA256=before,CFDRuns=0,physicalQualification=False)
    save(out/'report.json',result);print(json.dumps(result));return result['passed']

def check():
    from runOpenFoamProcess import run
    out=BASE/'layer_check01';out.mkdir(exist_ok=False);case=BASE/'layer_mesh01'
    assert json.loads((BASE/'layer_export01/report.json').read_text())['passed']
    r=json.loads((case/'report.json').read_text())
    assert all(sha(case/'constant/polyMesh'/n)==h for n,h in r['meshSHA256'].items())
    files=[p for d in ['constant','system'] for p in (case/d).rglob('*') if p.is_file()]
    files += [Path(__file__).resolve(),ROOT/'docs/validation/v208/LAYER_IMPORT_PLAN.md',case/'report.json']
    save(out/'admission.json',dict(seconds=60,parentSeconds=75,memoryMiB=2048,outputMiB=256,automaticRetry=0,inputSHA256={str(p):sha(p) for p in files}))
    shutil.copyfile(__file__,out/'executed_entry.py');start=time.monotonic()
    process=run(['checkMesh','-allGeometry','-allTopology'],case,out/'native',60,files+[out/'admission.json'],256,2048)
    text=(out/'native/stdout.txt').read_text()
    result=dict(passed=process['passed'] and 'Mesh OK.' in text,process=process,elapsedSeconds=time.monotonic()-start,
        originalQualityGatesUnchanged=True,CFDRuns=0,fullRegionalMeshAccepted=False,physicalQualification=False)
    save(out/'report.json',result);print(json.dumps(dict(passed=result['passed'],seconds=result['elapsedSeconds'],exit=process['nativeExitCode'])))
    return result['passed']

if __name__=='__main__':
    raise ValueError('v208原层导入和唯一原生质量检查均已结束；入口关闭，无CFD准入。')
