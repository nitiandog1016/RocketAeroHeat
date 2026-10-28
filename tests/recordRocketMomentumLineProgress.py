# -*- coding: utf-8 -*-
"""登记法向列的真实数值证据与报告恢复，拒绝把有界内部试算当作已解区域。"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,re,subprocess
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=ROOT/'docs/validation/v83/progress_momentum_lines.json';assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for b in iter(lambda:stream.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def identity(path):
 path=Path(path);return dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size)
def equal_numeric(a,b):
 if isinstance(a,dict):
  assert a.keys()==b.keys()
  for key in a:equal_numeric(a[key],b[key])
 else:np.testing.assert_array_equal(a,b)
names=['momentumlineunit01','momentumlineunit02','momentumlineunit03','momentumlineunit04','momentumlineunit05','momentumlinereal01','momentumlinereal02','momentumlinereal03','momentumlinereport01','momentumlinereport02','momentumlineui01'];processes={}
for name in names:
 file=VAL/(name+'_process.json');p=read(file);expected=4294967295 if name=='momentumlinereal03' else 1 if name=='momentumlinereport01' else 0
 assert p['status']=='terminal' and p['exitCode']==expected and all(p['sourceStableAfter'].values())
 processes[name]=dict(exitCode=p['exitCode'],sourceStable=True,engineProcessIDs=p['engineProcessIDs'],receipt=identity(file))
small=[read(VAL/(name+'_independent_momentum_lines.json')) for name in names[:5]];assert all(r['status'].startswith('passed-') for r in small)
real=read(VAL/'momentumlinereal03_independent_momentum_lines.json');assert real['status'].startswith('passed-') and real['innerMomentum']['completed'] and not real['innerMomentum']['internalTargetMet']
summary_file=VAL/'momentumlinereport02/momentum_lines_real_summary.json';assert summary_file.stat().st_size<128*1024;summary=read(summary_file)
assert not summary['wholeFiveFieldPreconditionerQualified'] and summary['driverReportingInterrupted'] and summary['summaryRecoveryNewAApplications']==summary['summaryRecoveryNewPreconditionerApplications']==0
stop=read(VAL/'momentumlinereal03/agent_reporting_stop_record.json');numeric=VAL/'momentumlinereal03/momentum_line_krylov_evidence.mat';assert sha(numeric)==stop['savedNumericalEvidenceSHA256']
archive=read(VAL/'momentumlinereal03/reporting_artifact_receipt.json');assert sha(archive['archivedPath'])==archive['sha256'];assert not (VAL/'momentumlinereal03/momentum_lines_real_summary.json').exists()
compact_unit=loadmat(VAL/'momentumlinereportunit01/momentum_summary_unit_evidence.mat',simplify_cells=True)['report'];assert compact_unit['status'].startswith('passed-') and max(compact_unit['summaryBytes'])<128*1024

# 预算修复仅涉及数组/元数据上界，完整缓存数值与修复前测试逐值相同。
before=loadmat(VAL/'momentumlineunit04/momentum_lines_unit_evidence.mat',simplify_cells=True);after=loadmat(VAL/'momentumlineunit05/momentum_lines_unit_evidence.mat',simplify_cells=True)
for a,b in zip(before['cases'],after['cases']):
 for key in ('inverseOrder','inverseLine','inverseNormInf','inverseCertification'):equal_numeric(a['factor'][key],b['factor'][key])
for key in ('frozenFactor','fallbackFactor'):
 for field in ('inverseOrder','inverseLine','inverseNormInf','inverseCertification'):equal_numeric(before[key][field],after[key][field])
factor_old=ROOT/'outputs/v83_validation/source_epochs/pre_momentum_metadata_budget_factorRocketMomentumLines.m';factor_current=ROOT/'internalFlow/cutCell/factorRocketMomentumLines.m'
assert sha(factor_old)==summary['sourceIdentity']['factorRocketMomentumLines']
old_text=factor_old.read_text(encoding='utf8');expected=old_text.replace('8*N*N*sourceFaces+8*n+sourceFaces*32+8192','8*N*N*sourceFaces+4*n+sourceFaces*64+32768').replace('assert(storage.bytes<=options.maximumFactorBytes,','assert(storage.bytes<=factorBound&&storage.bytes<=options.maximumFactorBytes,').replace('实际完整法向列因子超预算。','实际完整法向列因子超过登记上界或预算。')
assert expected==factor_current.read_text(encoding='utf8')
actual_bytes=int(re.search(r'MOMENTUM_LINES faces=19752 lineUnknowns=33 inverseDefect=\S+ bytes=(\d+)',(VAL/'momentumlinereal03_stdout.log').read_text(encoding='utf8',errors='replace')).group(1));new_bound=8*33*33*19752+4*217272+64*19752+32768;assert actual_bytes<=new_bound
audit_epoch=VAL/'source_epochs/passed_momentumlinereal03_audit_before_size_guard.py';assert sha(audit_epoch)==real['scriptSHA256']
current_audit=(ROOT/'tests/auditRocketMomentumLinesReal.py').read_bytes();without_guard=b''.join(line for line in current_audit.splitlines(keepends=True) if b'Use the explicit compact recovery report; never parse oversized basis JSON.' not in line)
assert hashlib.sha256(without_guard).hexdigest()==real['scriptSHA256']

view=VAL/'workbench22';data=read(view/'wall-velocity-study-data.json');research=read(VAL/'block_linear_research_summary_v6.json');assert data['linearResearch']==research and len(research['studies'])==9
row=research['studies'][-1];assert row['displayState']=='有限但内部未收敛' and row['innerRelativeCResidual'] is None and row['CApplications']==0 and not row['candidateUsable']
node=read(view/'node_wall_velocity_summary.json');assert node['status'].startswith('passed-') and node['CSVrows']==592560 and node['sectionSamples']==6624 and node['sameHostEvents']==15 and not node['browserVisualOrWebGLRuntimeVerified']
payload=[]
for case in data['cases']:
 rel=Path(case['resultsJSONURL'])
 for path in (rel,rel.parent/'case-data.js'):
  assert sha(view/path)==sha(VAL/'workbench21'/path);payload.append(identity(view/path))
assert len(payload)==6
ui=dict(status='passed-independent-nine-source-bound-research-rows-original-A-scope-and-unchanged-facet-data',processID=os.getpid(),workbench='outputs/v83_validation/workbench22/index.html',
 originalAResidualNotMislabelledAsCResidual=True,productionQualificationNotClaimed=True,sourceFacetPayload=payload,previousFullMATFieldAuditReusedByByteIdentity=True,browserVisualOrWebGLRuntimeVerified=False)
(view/'independent_momentum_line_study_view.json').write_text(json.dumps(ui,indent=2)+'\n',encoding='utf8')
old=read(ROOT/'docs/validation/v83/progress_consistent_schur.json');production=old['production']
for key in ('latestEndpoint','originalNSSource'):assert sha(ROOT/production[key]['path'])==production[key]['sha256']
assert subprocess.check_output(['git','rev-parse','--short','HEAD'],cwd=ROOT,text=True).strip()=='42601c9'
paths=[ROOT/'internalFlow/cutCell'/name for name in ('factorRocketMomentumLines.m','solveRocketMomentumLines.m','applyRocketMomentumGraph.m','solveRocketMomentumKrylov.m','rocketCompensatedInverseResidual.m','compactRocketMomentumKrylovReport.m')]
paths += [ROOT/'tests'/name for name in ('testRocketMomentumLines.m','studyRocketMomentumLines.m','testRocketMomentumSummary.m','recoverRocketMomentumStudySummary.m','diagnoseRocketMomentumLinePrecision.py','auditRocketMomentumLinesUnit.py','auditRocketMomentumLinesReal.py','exportRocketMomentumLineDiagnostics.py')]
out=dict(status='tested-certified-normal-or-original-CV-diagonal-preconditioner-with-unresolved-real-momentum-residual',version=83,recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid(),processes=processes,
 smallDomainAudits=small,realAudit=real,compactNumericalSummary=summary,compactSummaryReceipt=identity(summary_file),
 strictRejections=[dict(case=name,errorIdentifier=read(VAL/name/'momentum_lines_real_summary.json')['errorIdentifier'],factorAvailable=False) for name in ('momentumlinereal01','momentumlinereal02')],
 precisionReference=read(VAL/'momentumline_precision01/precision_reference.json'),reportingStop=stop,oversizedReportArchive=archive,
 summaryRecovery=dict(failedFirstSourcePathAttemptRetained=True,finalRecoveredWithoutNumericalRecomputation=True,compactUnitPassedInsideFailedFirstRecoveryBatch=True,compactSummaryBytes=summary_file.stat().st_size),
 metadataBudgetRepair=dict(numericalCachesExactlyUnchanged=True,onlyBudgetFormulaAndGuardChanged=True,actualRealFactorBytes=actual_bytes,newConservativeFactorBoundBytes=new_bound,
  originalNumericFactorSourceEpoch=identity(factor_old),passedRealAuditSourceEpoch=identity(audit_epoch),currentAuditAddsOnlyReportSizeGuard=True),
 codeSourceIdentity=[identity(p) for p in paths],userInterface=ui,nodeDataVerification=node,production=production,delivery=old['delivery'],
 workAccounting=dict(numericalRunFullAApplications=37,numericalRunLineInverseApplications=35,independentAuditFullAApplications=36,independentAuditLineMapApplications=35,
  summaryRecoveryNewAApplications=0,newSelectedJApplications=0,newActualNSEvaluations=0,newCoarseImageConstructions=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0),
 nextAction='Reuse the audited factor and saved real A/32 directions; do not reassemble normal lines, repeat global momentum ILU, restart stopped ILUTP, blindly extend iterations or scan fill/tau. Original A row-norm span is about 6.27e5 and current row-scaled residual is 1.0014, so first test an equivalent row-scale/preconditioner mapping on dense references, then one bounded same-RHS comparison retaining original-unit A residual acceptance. Scaling is unproven as a convergence repair. If still weak, assess cross-line coupling and the field grouping; a bounded internal solve is not a fixed Schur inverse. Preserve all CVs, original C/NS and all direction/temperature/descent/convergence gates, regional verification and scientific UI.',
 goalStatus='active',wholeGoalComplete=False,recorderSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],manifest=identity(target),researchRows=9)))
