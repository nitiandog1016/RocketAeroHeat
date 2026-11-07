"""从既有正确原生查询回执收口物性，检查一次修正后的真三维网格。"""
from pathlib import Path
import json,time,shutil
import numpy as np
from prepareOpenFoamControl import ROOT,BASE,save,sha
from openFoamAirThermo import evaluate,contract,R,HF

def main():
    from runOpenFoamProcess import run
    start=time.monotonic();out=BASE/'prepared_audit01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    query=BASE/'native_finish01/query';q=json.loads((query/'report.json').read_text());assert q['passed']
    text=(query/'stdout.txt').read_text();rows=np.array([[float(x) for x in s.split()[1:]] for s in text.splitlines() if s.startswith('ROW ')])
    meta=next(s.split() for s in text.splitlines() if s.startswith('R '))
    assert len(rows)==9 and abs(float(meta[1])/R-1)<2e-15 and abs(float(meta[3])-HF)<2e-9
    g=evaluate(5000,rows[:,0]);errors={}
    for i,key in enumerate(['cpJkgK','cvJkgK','enthalpyJkg','sensibleEnthalpyJkg','internalEnergyJkg','densityKgm3','soundSpeedMps']):
        errors[key]=float(np.max(abs(rows[:,i+1]-g[key])/np.maximum(1,abs(g[key]))));assert errors[key]<2e-10
    inv=float(np.max(abs(rows[:,-1]-rows[:,0])));assert inv<2e-7
    case=BASE/'control03';prep=json.loads((BASE/'preparation_repaired.json').read_text());assert prep['passed']
    for p,h in prep['inputSHA256'].items():assert sha(p)==h
    for rel in ['src/meshCheck/primitiveMeshCheck/primitiveMeshCheck.C','src/meshCheck/polyMeshCheck/polyMeshCheck.C']:
        p=Path('/opt/openfoam13')/rel;shutil.copyfile(p,out/p.name)
    native=run(['checkMesh','-allTopology','-allGeometry'],case,out/'mesh',60,
        [__file__,*[p for p in (case/'constant/polyMesh').iterdir()]],96)
    text=(out/'mesh/stdout.txt').read_text();assert native['passed'] and 'Mesh OK.' in text
    report=dict(passed=True,queryReceiptSHA256=sha(query/'report.json'),queryLogSHA256=sha(query/'stdout.txt'),
        propertyMaxRelativeErrors=errors,inverseTemperatureMaxAbsK=inv,contract=contract(),
        geometry=prep['geometry'],mesh=native,originalMeshCheckOverallPassed=False,
        correctedMeshCheckOverallPassed=True,CFDRuns=0,engineBuilds=0,
        elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','elapsedSeconds','geometry','propertyMaxRelativeErrors']}))

if __name__=='__main__':main()
