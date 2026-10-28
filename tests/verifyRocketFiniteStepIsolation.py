"""核对独立工程重算与冻结六步证据，不把完整Q24回放当成整箭多步重算。"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v84_validation'
STAGE=OUT/'delivery/isolated/RocketAeroHeat'
RUN=STAGE/'outputs/isolated_validation'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()

def main():
    proof=json.loads((RUN/'isolated_package.json').read_text('utf8'))
    assert proof['status'].startswith('passed-isolated-package')
    Path(proof['resolvedNSFile']).relative_to(STAGE)
    assert proof['completeQ24Replay']['rateAndLocalInfoBitwiseIdentical']
    assert not proof['fullMultistepComparisonRecomputedInIsolation']
    cases=[]
    for original,recomputed in [('patch02','comparison'),('patch03','matched')]:
        a=json.loads((OUT/original/'patch_comparison.json').read_text('utf8'))
        b=json.loads((RUN/recomputed/'patch_comparison.json').read_text('utf8'))
        frozen=loadmat(OUT/original/'patch_comparison.mat',simplify_cells=True)
        fresh=loadmat(RUN/recomputed/'patch_comparison.mat',simplify_cells=True)
        for branch in ['baseline','guarded']:
            assert a[branch]['acceptedSteps']==b[branch]['acceptedSteps']==6
            assert a[branch]['residualEvaluations']==b[branch]['residualEvaluations']
            for x,y in zip(a[branch]['history'],b[branch]['history']):
                assert {k:v for k,v in x.items() if k!='elapsedSeconds'}=={k:v for k,v in y.items() if k!='elapsedSeconds'}
            for field in ['Q','rate']:
                np.testing.assert_array_equal(frozen[branch][field],fresh[branch][field])
            for x,y in zip(frozen[branch]['directions'],fresh[branch]['directions']):
                for field in ['Q','direction','jv']:
                    np.testing.assert_array_equal(x[field],y[field])
            cases.append({'case':original,'branch':branch,'sixStepHistoryExactlyReplayed':True,
                          'allSavedStatesDirectionsAndRatesBitwiseIdentical':True,'residualEvaluations':b[branch]['residualEvaluations']})
    report={'status':'passed-isolated-six-step-state-direction-rate-and-complete-Q24-replay-verification',
            'cases':cases,'completeQ24ReplaySeconds':proof['completeQ24ReplaySeconds'],
            'fullMultistepComparisonRecomputedInIsolation':False,'sourceSHA256':sha(Path(__file__)),
            'isolatedProofSHA256':sha(RUN/'isolated_package.json'),'newProductionSolvedAreaM2':0}
    (OUT/'delivery/isolated_comparison.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    (RUN/'independent_comparison.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report))

if __name__=='__main__':main()
