# -*- coding: utf-8 -*-
"""核对像列归一研究的F/父C资格分离、旧18行及全部原面字段字节保持。"""
from pathlib import Path
import json
import os
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


folder=VAL/'workbench29'
target=folder/'independent_full_F_image_scale_study_view.json'
assert not target.exists()
d=read(folder/'wall-velocity-study-data.json');old=read(VAL/'workbench28/wall-velocity-study-data.json')
r,p=d.pop('linearResearch'),old.pop('linearResearch')
assert d==old and r['studies'][:18]==p['studies'] and len(r['studies'])==19
assert r==read(VAL/'block_linear_research_summary_v13.json')
bindings={}
for row in r['studies']:
    file=Path(row['evidenceFile']);file=file if file.is_absolute() else ROOT/file
    assert sha(file)==row['evidenceSHA256'] and not row['candidateUsable']
    assert row['newActualNSEvaluations']==row['newProductionSolvedAreaM2']==0
    bindings[file.relative_to(ROOT).as_posix()]=sha(file)
row=r['studies'][-1];meta=r['fullFImageScaleDiagnostics']
assert meta['caseID']==row['caseID'] and meta['algorithm']=='LSQR'
assert row['innerRelativeCResidual'] is None and row['basisOrthogonalityError'] is None and row['modelRelationRelativeDifference'] is None
assert row['arnoldiSteps']==0 and row['CApplications']==1
assert meta['fullFInternalTargetMet']==(meta['finalFullFResidual']<=1e-3)
assert meta['parentOriginalCTargetMet']==(meta['parentOriginalCResidual']<=1e-3)
assert meta['allRealColumnsIndependentlyCompared']==187749*5
assert sha(meta['auditFile'])==meta['auditSHA256']
assert not r['physicalStateChanged'] and not r['productionQualified'] and not meta['originalNSOrRegionalQualificationProven']
payload=[]
for case in d['cases']:
    for name in ('results.json','case-data.js'):
        relative=Path(case['resultsJSONURL']).parent/name
        before,after=VAL/'workbench28'/relative,folder/relative
        assert sha(before)==sha(after)
        payload.append(dict(path=after.relative_to(ROOT).as_posix(),sha256=sha(after),bytes=after.stat().st_size))
assert len(payload)==6
node=read(folder/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed-') and node['CSVrows']==592560 and node['sectionSamples']==6624 and node['sameHostEvents']==15
assert not node['browserVisualOrWebGLRuntimeVerified']
binding=loadmat(VAL/'fullFimagescaleui01/linear_study_binding_evidence.mat',simplify_cells=True)['report'];assert all(binding['rejected'])
process=read(VAL/'fullFimagescaleui01_process.json')
assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
out=dict(status='passed-independent-nineteen-source-bound-rows-image-scale-F-parent-qualification-and-unchanged-facet-data',
    processID=os.getpid(),researchRows=19,previousEighteenRowsExactlyPreserved=True,
    originalTenFieldsAndNSHistoryUnchanged=True,fullFAndParentAndNSQualificationKeptSeparate=True,
    originalFNormAndRightHandSideUnchanged=True,LSQRArnoldiMetricsNotApplicable=True,
    allThreeFalseSourceOrQualificationExportsRejected=True,sourceFacetPayload=payload,
    previousFullMATFieldAuditReusedByByteIdentity=True,browserVisualOrWebGLRuntimeVerified=False,
    newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','sourceFacetPayload')}))
