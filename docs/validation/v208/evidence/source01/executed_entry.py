"""固定来源的出口与热物性接口盘点；只复制小源码，不构建或推进流场。"""
from pathlib import Path
import json,shutil,time
from prepareOpenFoamControl import ROOT,sha,save

def main():
    start=time.monotonic();out=ROOT/'outputs/v208_validation/source01';out.mkdir(parents=True,exist_ok=False)
    engine=Path('/opt/openfoam13');base=engine/'src/finiteVolume/fields/fvPatchFields/derived'
    names=[p.name for p in base.iterdir() if any(k in p.name.lower() for k in ['wave','reflect','acoustic','advect','character'])]
    files=[]
    for name in ['advective','waveTransmissive','fixedMeanOutletInlet','inletOutlet','zeroGradient']:
        directory=base/name
        if directory.exists():files += [p for p in directory.iterdir() if p.suffix in ['.H','.C']]
    for sub in ['src/thermophysicalModels/basic/basicThermo','src/thermophysicalModels/basic/fluidThermo',
                'src/thermophysicalModels/basic/heThermo','src/thermophysicalModels/basic/psiThermo',
                'src/thermophysicalModels/basic/hePsiThermo']:
        directory=engine/sub
        if directory.exists():files += [p for p in directory.iterdir() if p.suffix in ['.H','.C']]
    records=[]
    for p in files:
        if time.monotonic()-start>60:raise TimeoutError('source budget')
        rel=p.relative_to(engine);dest=out/'source'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        records.append(dict(path=str(p),relativePath=rel.as_posix(),sha256=sha(p),bytes=p.stat().st_size))
    shutil.copyfile(__file__,out/'executed_entry.py')
    save(out/'report.json',dict(passed=True,sources=records,relatedNativeBoundaryNames=names,
        fixedEngineCommit='18870c24d21c6b982e2cdec27b2f59738cca5f90',elapsedSeconds=time.monotonic()-start,
        CFDRuns=0,builds=0,physicalQualification=False))
    print(json.dumps(dict(files=len(records),bytes=sum(r['bytes'] for r in records),names=names,seconds=time.monotonic()-start)))

if __name__=='__main__':main()
