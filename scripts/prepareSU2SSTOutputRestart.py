"""以新格式原生检查点做两步小型兼容验证；不重做平板或完整域物理解。"""
import os,json
from prepareSU2CheckpointContinuation import ROOT,read
from prepareSU2SSTCase import set_value
from prepareSU2SSTCheckpoint import prepare

def main():
    source=ROOT/'outputs/v126_validation/format01';folder=ROOT/'outputs/v126_validation/restart02'
    assert not folder.exists();folder.mkdir(parents=True)
    os.link(str(source/'mesh.su2'),str(folder/'mesh.su2'))
    text=set_value((source/'case.cfg').read_text(encoding='utf-8'),'ITER',2)
    prepare(folder,source,text,read(source/'sst_contract.json'),read(source/'case_manifest.json'))
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=30,
        dryRun=False,sourceFiles=['scripts/prepareSU2SSTOutputRestart.py','docs/validation/v126/PLAN.md'])
    (folder/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(directory=str(folder),steps=2)))

if __name__=='__main__':main()
