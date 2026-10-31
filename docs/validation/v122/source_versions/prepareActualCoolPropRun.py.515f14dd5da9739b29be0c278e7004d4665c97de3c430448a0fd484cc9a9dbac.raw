"""复用完整真实网格和原投影，准备唯一20步变比热任务；不重做网格或直接复用旧能量。"""
from pathlib import Path
import hashlib,json,os,re,sys,uuid
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CoolPropCase import prepare
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=ROOT/'outputs/v120_validation/study01/control';folder=ROOT/'outputs/v122_validation/actual01'
    assert not folder.exists();folder.mkdir(parents=True)
    for name in ['mesh.su2','projection.mat']:os.link(str(source/name),str(folder/name))
    cfg=(source/'case.cfg').read_text()
    for key,value in [('RESTART_SOL','NO'),('ITER','20'),('OUTPUT_WRT_FREQ','10')]:
        cfg,n=re.subn(r'^'+key+r'=.*$',key+'= '+value,cfg,flags=re.M);assert n==1
    for key in ['SOLUTION_FILENAME','READ_BINARY_RESTART']:cfg,n=re.subn(r'^'+key+r'=.*\n','',cfg,flags=re.M);assert n==1
    (folder/'case.cfg').write_text(cfg)
    manifest=json.loads((source/'case_manifest.json').read_text());manifest.pop('continuation',None);manifest.pop('initialization',None)
    manifest['flight']['iterations']=20;manifest['flight']['outputEvery']=10
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    r=prepare(folder,source)
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=600,dryRun=False,jobToken=str(uuid.uuid4()),
        sourceFiles=[str(p.relative_to(ROOT)) for p in [Path(__file__).resolve(),ROOT/'scripts/prepareSU2CoolPropCase.py',
            ROOT/'docs/validation/v122/PLAN.md',folder/'projection.mat',folder/'case_manifest.json',folder/'farfield_audit.json']])
    (folder/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['nodeCount','errors','sourceCumulativeIteration','EOSMach']},indent=2))
if __name__=='__main__':main()
