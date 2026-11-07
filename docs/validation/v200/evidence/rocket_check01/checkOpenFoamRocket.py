"""原生完整几何检查和边界源绑定；不运行真实箭体流场。"""
from pathlib import Path
import json,shutil,time
from prepareOpenFoamControl import ROOT,sha,save
from runOpenFoamProcess import run


def main():
    start=time.monotonic();base=ROOT/'outputs/v200_validation';case=base/'rocket02';out=base/'rocket_check01';out.mkdir(exist_ok=False)
    prep_path=base/'rocket_preparation02/report.json';prep=json.loads(prep_path.read_text());assert prep['passed'] and not prep['executionAuthorized']
    for p,h in prep['inputSHA256'].items():assert sha(str(p).replace('E:\\','/mnt/e/').replace('\\','/'))==h
    for name in ['checkOpenFoamRocket.py','runOpenFoamProcess.py']:shutil.copyfile(ROOT/'scripts'/name,out/name)
    # 模块一：实际freestream继承inletOutlet；用通量符号选择来流或零梯度。
    src=Path('/opt/openfoam13/src/finiteVolume/fields/fvPatchFields')
    sources=[]
    for name in ['freestreamFvPatchField.C','inletOutletFvPatchField.C']:
        hits=[p for p in src.rglob(name) if 'lnInclude' not in p.parts];assert len(hits)==1
        p=hits[0];shutil.copyfile(p,out/name);sources.append(p)
    a=sources[0].read_text();b=sources[1].read_text()
    assert 'inletOutletFvPatchField' in a and 'freestreamValue' in a
    assert '1.0 - pos0(phip)' in b
    locks=[Path(__file__),prep_path,ROOT/'docs/validation/v200/ROCKET_ENTRY_PLAN.md',*sources,*[p for p in case.rglob('*') if p.is_file()]]
    process=run(['checkMesh','-allTopology','-allGeometry'],case,out/'native',60,locks,96,2048)
    text=(out/'native/stdout.txt').read_text();passed=process['passed'] and 'Mesh OK.' in text and 'Failed' not in text
    result=dict(passed=passed,process=process,geometry=prep['geometry'],originalWallExact=True,originalWallFaces=19752,
        originalWallAreaM2=prep['originalWallAreaM2'],boundarySourceSHA256={str(p):sha(p) for p in sources},
        originalRocketCFDRuns=0,engineBuilds=0,executionAuthorized=False,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(dict(passed=passed,seconds=result['elapsedSeconds'],stop=process['stopReason'])))


if __name__=='__main__':main()
