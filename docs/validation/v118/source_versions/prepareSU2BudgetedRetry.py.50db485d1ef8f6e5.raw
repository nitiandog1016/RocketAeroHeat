"""实际预算失败后显式重试；复用原初态/网格，改变迭代预算及具名保存间隔。"""
from pathlib import Path
import argparse,json,hashlib,os,re,uuid
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--target',type=Path,required=True)
    ap.add_argument('--iterations',type=int,required=True);ap.add_argument('--every',type=int,required=True);a=ap.parse_args()
    src=a.source.resolve();dst=a.target.resolve();src.relative_to(ROOT/'outputs');dst.relative_to(ROOT/'outputs')
    r=json.loads((src/'engine_process.json').read_text(encoding='utf-8'))
    assert r['status']=='terminal' and r['cancelled'] and all(r['inputsStable'].values())
    assert 1<=a.every<a.iterations and a.iterations<600 and not dst.exists()
    projected=r['elapsedSeconds']/max(1,r['progress']['Inner_Iter']+1)*a.iterations
    assert projected<.8*720,'实测预算不足以安全完成，不能盲重试'
    for name in ['mesh.su2','initial_flow.csv','projection.mat']:
        assert sha(src/name)==r['inputSHA256'][str(src/name)]
    dst.mkdir()
    for name in ['mesh.su2','initial_flow.csv','projection.mat','initialization_manifest.json']:os.link(str(src/name),str(dst/name))
    cfg=(src/'case.cfg').read_text(encoding='utf-8');before=cfg
    for key,value in [('ITER',a.iterations),('OUTPUT_WRT_FREQ',a.every)]:
        cfg,n=re.subn(r'^'+key+r'=.*$',key+'= '+str(value),cfg,flags=re.M);assert n==1
    (dst/'case.cfg').write_text(cfg,encoding='utf-8')
    m=json.loads((src/'case_manifest.json').read_text(encoding='utf-8'));m['flight'].update(iterations=a.iterations,outputEvery=a.every)
    m['budgetRetry']=dict(originalCancelledCase=str(src),originalReceiptSHA256=sha(src/'engine_process.json'),
        restartFromOriginalInitialGuess=True,notRecoveredFromCancelledIteration=True,
        projectedSeconds=projected,physicalInputsUnchanged=True,physicsQualified=False)
    (dst/'case_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    sources=[str(dst/'projection.mat'),str(dst/'case_manifest.json'),str(src/'engine_process.json'),str(src/'case.cfg'),
        'scripts/prepareSU2BudgetedRetry.py','writeSU2HybridInitialState.m','outputs/v118_validation/budget_revision.json']
    request=dict(workDirectory=str(dst),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=720,
        dryRun=False,jobToken=str(uuid.uuid4()),sourceFiles=sources)
    (dst/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n')
    print(json.dumps(dict(projectedSeconds=projected,iterations=a.iterations,every=a.every,request=str(dst/'engine_request.json'))))
if __name__=='__main__':main()
