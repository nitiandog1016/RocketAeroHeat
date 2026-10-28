# -*- coding: utf-8 -*-
"""核对新增研究摘要的原证据绑定，及旧/新工作台原面场量和历史的字节一致性。"""
from pathlib import Path
import hashlib,json,os,sys
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';folder=Path(sys.argv[1]);folder=folder if folder.is_absolute() else ROOT/folder
target=folder/'independent_block_study_view.json';assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
summaryFile=Path(sys.argv[2]) if len(sys.argv)>2 else VAL/'block_linear_research_summary.json';summaryFile=summaryFile if summaryFile.is_absolute() else ROOT/summaryFile
data=read(folder/'wall-velocity-study-data.json');old=read(VAL/'workbench16/wall-velocity-study-data.json');research=data.pop('linearResearch');assert data==old
assert research==read(summaryFile) and research['sourceCase']=='krylov17' and research['sourceEvidenceSHA256']==sha(VAL/'krylov17/wall_velocity_krylov_evidence.mat')
assert not research['physicalStateChanged'] and not research['productionQualified'] and research['newProductionSolvedAreaM2']==0
bindings={};unchanged=[]
for c in data['cases']:
    for name in ('results.json','case-data.js'):
        relative=Path(c['resultsJSONURL']).parent/name;before=VAL/'workbench16'/relative;after=folder/relative
        assert sha(before)==sha(after);unchanged.append(str(relative));bindings[str(after.relative_to(ROOT))]=sha(after)
for s in research['studies']:
    evidence=Path(s['evidenceFile']);evidence=evidence if evidence.is_absolute() else ROOT/evidence;assert sha(evidence)==s['evidenceSHA256']
    assert not s['candidateUsable'] and s['newActualNSEvaluations']==0 and s['newProductionSolvedAreaM2']==0
    if s['modelRelationRelativeDifference'] is not None:assert s['modelRelationRelativeDifference']>=s['modelRelationThreshold']==1e-8
    else:assert s['displayState'] in ('分解拒绝','单次解不稳定','原分解拒绝（已定位）','零主元拒绝') and s['basisOrthogonalityError'] is None
    bindings[str(evidence.relative_to(ROOT))]=s['evidenceSHA256']
if 'failedCellDetails' in research:
    defect=research['failedCellDetails']['localInverseDefectValue'];assert research['studies'][-1]['innerRelativeCResidual'] is None
    if defect is not None:
        receipt=Path(research['pivotFailureCheckpoint']['receiptFile']);receipt=receipt if receipt.is_absolute() else ROOT/receipt
        assert sha(receipt)==research['pivotFailureCheckpoint']['receiptSHA256'] and read(receipt)['failureSummary']['maximumInverseDefect']==defect and defect>=1e-8
node=read(folder/'node_wall_velocity_summary.json');assert node['sameHostEvents']==15 and node['CSVrows']==592560 and node['sectionSamples']==6624 and not node['browserVisualOrWebGLRuntimeVerified']
for file in [folder/'wall-velocity-study-data.json',folder/'index.html',folder/'rocket-wall-velocity.js',folder/'node_wall_velocity_summary.json',VAL/'workbench16/independent_wall_velocity_view.json',summaryFile]:bindings[str(file.relative_to(ROOT))]=sha(file)
out=dict(status='passed-independent-research-source-binding-and-unchanged-wall-fields-history-identity',processID=os.getpid(),
         unchangedSourceFacetPayloadFiles=unchanged,separateRejectedInternalStudyRows=len(research['studies']),originalNSAndAcceptanceHistoryExactlyUnchanged=True,
         priorIndependentFullMATFieldAudit='outputs/v83_validation/workbench16/independent_wall_velocity_view.json',priorFieldAuditReusedByExactPayloadIdentity=True,
         newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,browserVisualOrWebGLRuntimeVerified=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
