# -*- coding: utf-8 -*-
"""核对新完整父C行、各级资格分离以及旧十九行/原面结果的字节保持。"""
from pathlib import Path
import json,os
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha


def read(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))


folder=VAL/'workbench30';target=folder/'independent_flexible_reduction_study_view.json';assert not target.exists()
d=read(folder/'wall-velocity-study-data.json');old=read(VAL/'workbench29/wall-velocity-study-data.json')
r,p=d.pop('linearResearch'),old.pop('linearResearch')
assert d==old and r['studies'][:19]==p['studies'] and len(r['studies'])==20
assert r==read(VAL/'block_linear_research_summary_v14.json')
meta=r['flexibleReductionDiagnostics'];bindings={}
for row in r['studies']:
    path=Path(row['evidenceFile']);path=path if path.is_absolute() else ROOT/path
    assert sha(path)==row['evidenceSHA256']
    assert row['newActualNSEvaluations']==row['newProductionSolvedAreaM2']==0
    expected=bool(meta['internalCTargetMet']) if row['caseID']==meta['caseID'] else False
    assert row['candidateUsable']==expected
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
row=r['studies'][-1]
assert row['caseID']==meta['caseID'] and row['innerRelativeCResidual']==meta['relativeOriginalCResidual']
assert meta['algorithm']=='FGMRES' and meta['internalCTargetMet']==(meta['relativeOriginalCResidual']<=1e-3)
assert row['basisOrthogonalityError']<1e-8 and row['modelRelationRelativeDifference']<1e-8
assert meta['allRealOuterArnoldiColumnsIndependentlyVerified'] and not meta['independentAllRealInnerLSQRColumnsClaimed']
assert sha(meta['auditFile'])==meta['auditSHA256'] and not meta['oldSchurIdentityClaimed']
assert not r['physicalStateChanged'] and not r['productionQualified'] and not meta['originalNSOrRegionalQualificationProven']
payload=[]
for case in d['cases']:
    for name in ('results.json','case-data.js'):
        relative=Path(case['resultsJSONURL']).parent/name
        before,after=VAL/'workbench29'/relative,folder/relative
        assert sha(before)==sha(after)
        payload.append(dict(path=after.relative_to(ROOT).as_posix(),sha256=sha(after),bytes=after.stat().st_size))
assert len(payload)==6
node=read(folder/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows']==592560 and node['sectionSamples']==6624 and node['sameHostEvents']==15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding=loadmat(VAL/'flexiblereductionui01/linear_study_binding_evidence.mat',simplify_cells=True)['report'];assert all(binding['rejected'])
process=read(VAL/'flexiblereductionui01_process.json')
assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
out=dict(status='passed-independent-twenty-source-bound-rows-flexible-C-qualification-and-unchanged-facet-data',processID=os.getpid(),researchRows=20,
    previousNineteenRowsExactlyPreserved=True,originalTenFieldsAndNSHistoryUnchanged=True,
    originalCInternalQualificationNotPromotedToNSOrRegion=True,allRealOuterColumnAuditBound=True,
    allThreeFalseSourceOrQualificationExportsRejected=True,sourceFacetPayload=payload,
    previousFullMATFieldAuditReusedByByteIdentity=True,browserVisualOrWebGLRuntimeVerified=False,
    newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','sourceFacetPayload')}))
