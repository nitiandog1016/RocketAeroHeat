"""复用完整SST终态的一次空间诊断；数值/物理保持，源149与新局部0分列。"""
import copy,json,os
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCheckpoint import prepare as checkpoint
from prepareSU2SSTCase import set_value
from su2SSTInterface import SPATIAL_OUTPUT


def main():
    source=ROOT/'outputs/v127_validation/control';folder=ROOT/'outputs/v128_validation/spatial01'
    assert not folder.exists();folder.mkdir(parents=True)
    plan=ROOT/'docs/validation/v128/PLAN.md';assert plan.is_file()
    text=(source/'case.cfg').read_text(encoding='utf-8-sig')
    for key,value in dict(ITER=1,OUTPUT_WRT_FREQ=1,SCREEN_WRT_FREQ_INNER=1,
        VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL, TIMESTEP)').items():
        text=set_value(text,key,value)
    os.link(str(source/'mesh.su2'),str(folder/'mesh.su2'))
    contract=copy.deepcopy(read(source/'sst_contract.json'));contract.pop('numericalStudy',None)
    contract.update(spatialOutputContract=SPATIAL_OUTPUT,spatialDiagnosticOnly=True)
    manifest=copy.deepcopy(read(source/'case_manifest.json'))
    manifest['sourceRunControlsFromConfig']={k:config(source/'case.cfg')[k] for k in ['ITER','OUTPUT_WRT_FREQ']}
    manifest['flight'].update(iterations=1,outputEvery=1,namedCheckpoints=True)
    checkpoint(folder,source,text,contract,manifest)
    init=read(folder/'sst_initialization.json')
    init.update(diagnosticOnly=True,diagnosticPurpose='spatial-implicit-RHS-and-pseudotime')
    (folder/'sst_initialization.json').write_text(json.dumps(init,indent=2)+'\n',encoding='utf-8')
    contract['initializationSHA256']=sha(folder/'sst_initialization.json')
    (folder/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    manifest.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'))
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,
        maximumSeconds=180,dryRun=False,sourceFiles=['scripts/prepareSU2SSTSpatialDiagnostic.py','docs/validation/v128/PLAN.md'])
    (folder/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(folder=str(folder),sourceLocalIteration=149,newLocalIteration=0,updates=1)))


if __name__=='__main__':main()
