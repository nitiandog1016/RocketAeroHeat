"""复用v116原24层网格/初态，建立有显式模式和具名检查点的新研究目录。"""
from pathlib import Path
import argparse,json,hashlib,os,uuid,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v117_validation'
SOURCE=ROOT/'outputs/v116_validation/study01/control'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['diagnostic','apply']);a=p.parse_args()
    if a.mode=='apply':
        gate=json.loads((OUT/'diagnostic_gate.json').read_text())
        assert gate['passed'] and gate['eligibleTargetCount']==4
        assert sha(OUT/'diagnostic/rah_frozen_directions.csv')==gate['directionsSHA256']
    receipt=json.loads((SOURCE/'engine_process.json').read_text())
    for name in ['mesh.su2','initial.csv','projection.mat']:
        path=SOURCE/name;assert sha(path)==receipt['inputSHA256'][str(path)]
    dest=OUT/a.mode;assert not dest.exists();dest.mkdir()
    for name in ['mesh.su2','initial.csv','projection.mat']:os.link(str(SOURCE/name),str(dest/name))
    steps=2 if a.mode=='diagnostic' else 600
    text=(SOURCE/'case.cfg').read_text()
    updates={'ITER':str(steps),'OUTPUT_WRT_FREQ':str(steps),'WRT_RESTART_OVERWRITE':'NO','WRT_SURFACE_OVERWRITE':'NO'}
    for key,value in updates.items():
        text,count=re.subn(r'^'+key+r'=.*$',key+'= '+value,text,flags=re.M);assert count==1
    (dest/'case.cfg').write_text(text,encoding='utf-8')
    manifest=json.loads((SOURCE/'case_manifest.json').read_text(encoding='utf-8'))
    manifest.update(requiredRuntimeFlavor='rah-frozen-energy-v1',frozenEnergyMode=a.mode,
        physicsQualified=False,study='v117 frozen-energy mechanism; original geometry and threshold retained')
    manifest['flight'].update(iterations=steps,outputEvery=steps,namedCheckpoints=True)
    (dest/'case_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    request=dict(workDirectory=str(dest),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=720,
        dryRun=False,jobToken=str(uuid.uuid4()),sourceFiles=[
            str(dest/'projection.mat'),str(dest/'case_manifest.json'),'scripts/prepareSU2FrozenEnergyStudy.py',
            'third_party/su2_thermal/RAHFrozenEnergyRecovery.hpp','third_party/su2_thermal/RAHFrozenEnergyDiagnostics.hpp',
            'third_party/su2_thermal/FrozenEnergyHook.cpp.inc'])
    if a.mode=='apply':request['sourceFiles'].append('outputs/v117_validation/diagnostic_gate.json')
    (dest/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n')
    print(str(dest/'engine_request.json'))
if __name__=='__main__':main()
