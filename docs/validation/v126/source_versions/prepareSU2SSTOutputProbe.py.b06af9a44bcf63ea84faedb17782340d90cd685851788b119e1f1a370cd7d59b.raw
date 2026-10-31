"""新格式的最小三步任务；复用旧平板输入但不修改旧场或重做物理验证。"""
from pathlib import Path
import os,json
from prepareSU2CheckpointContinuation import ROOT,sha,read
from prepareSU2SSTCase import set_value
from su2SSTInterface import PHASE_OUTPUT


def main():
    source=ROOT/'outputs/v124_validation/plate01/plate';target=ROOT/'outputs/v126_validation/format01'
    assert not target.exists();r=read(source/'engine_process.json')
    assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
    for name in ['case.cfg','mesh.su2','case_manifest.json','sst_contract.json']:
        assert sha(source/name)==r['inputSHA256'][str(source/name)]
    target.mkdir(parents=True);os.link(str(source/'mesh.su2'),str(target/'mesh.su2'))
    text=(source/'case.cfg').read_text(encoding='utf-8-sig')
    for k,v in dict(ITER=3,OUTPUT_WRT_FREQ=1,SCREEN_WRT_FREQ_INNER=1,WRT_RESTART_COMPACT='NO').items():text=set_value(text,k,v)
    (target/'case.cfg').write_text(text,encoding='utf-8')
    contract=read(source/'sst_contract.json');contract.update(outputStateContract=PHASE_OUTPUT,configurationSHA256=sha(target/'case.cfg'))
    (target/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    manifest=read(source/'case_manifest.json');manifest['sst']=contract;manifest['sstContractSHA256']=sha(target/'sst_contract.json')
    manifest['formatStudy']=dict(sourceDirectory=str(source),sourceReceiptSHA256=sha(source/'engine_process.json'),
        coldStart=True,updates=3,physicsQualified=False)
    (target/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    request=dict(workDirectory=str(target),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=30,dryRun=False,
        sourceFiles=['scripts/prepareSU2SSTOutputProbe.py','docs/validation/v126/PLAN.md'])
    (target/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(directory=str(target),steps=3)))

if __name__=='__main__':main()
