# -*- coding: utf-8 -*-
"""核对完整父C的第21行、原NS探针资格分离与原面字节身份。"""
from pathlib import Path
import json,os
from screenRocketFullFLeastSquares import ROOT,VAL,sha

def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))

folder=VAL/'workbench31';target=folder/'independent_parent_least_squares_view.json';assert not target.exists()
d=read(folder/'wall-velocity-study-data.json');old=read(VAL/'workbench30/wall-velocity-study-data.json')
r,p=d.pop('linearResearch'),old.pop('linearResearch')
assert d==old and len(r['studies'])==21 and r['studies'][:20]==p['studies']
assert r==read(VAL/'block_linear_research_summary_v15.json')
meta=r['parentLeastSquaresDiagnostics'];probe=r['sameQOriginalNSDirectionProbe'];row=r['studies'][-1]
assert row['caseID']==meta['caseID']=='parentlsqrreal01' and meta['algorithm']=='LSQR'
assert row['innerRelativeCResidual']==meta['relativeOriginalCResidual'] and row['basisOrthogonalityError'] is None
assert row['modelRelationRelativeDifference'] is None and not row['candidateUsable']
assert not meta['internalCTargetMet'] and meta['relativeOriginalCResidual']>meta['originalCTarget']==1e-3
assert probe['directionTargetMet'] and probe['actualTotalDirectionDefect']<probe['originalTarget']==.1
assert not probe['physicalStateChanged'] and not probe['regionalQualificationProven']
assert probe['newActualNSEvaluations']==2 and probe['newFlowAcceptances']==probe['newProductionSolvedAreaM2']==0
assert not meta['independentAllRealBidiagonalizationColumnsClaimed']
assert not meta['oldSchurUsed'] and not meta['fixedCRepresentativesUsed']
bindings={}
for item in r['studies']:
    path=Path(item['evidenceFile']);path=path if path.is_absolute() else ROOT/path
    assert sha(path)==item['evidenceSHA256'] and item['newActualNSEvaluations']==item['newProductionSolvedAreaM2']==0
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
for item in (meta,probe):
    path=Path(item['auditFile']);assert sha(path)==item['auditSHA256']
    assert read(path)['status'].startswith('passed')
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
payload=[]
for case in d['cases']:
    for name in ('results.json','case-data.js'):
        relative=Path(case['resultsJSONURL']).parent/name
        before,after=VAL/'workbench30'/relative,folder/relative
        assert sha(before)==sha(after)
        payload.append(dict(path=after.relative_to(ROOT).as_posix(),sha256=sha(after)))
node=read(folder/'node_wall_velocity_summary.json')
assert node['status'].startswith('passed') and node['CSVrows']==592560 and node['sectionSamples']==6624 and node['sameHostEvents']==15
assert not node['browserVisualOrWebGLRuntimeVerified']
process=read(VAL/'parentlsqrui01_process.json')
assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
out=dict(status='passed-independent-21-source-bound-rows-parent-C-and-NS-qualification-separation',processID=os.getpid(),
    researchRows=21,previousTwentyRowsExactlyPreserved=True,sourceFacetPayload=payload,
    originalTenFieldsAndNSHistoryUnchanged=True,previousFullMATFieldAuditReusedByByteIdentity=True,
    internalCTargetStillNotMet=True,separateNSDirectionGatePassed=True,newProductionSolvedAreaM2=0,
    browserVisualOrWebGLRuntimeVerified=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','sourceFacetPayload')}))
