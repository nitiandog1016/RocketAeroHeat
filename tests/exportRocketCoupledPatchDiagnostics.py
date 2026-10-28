# -*- coding: utf-8 -*-
"""接入等价行尺度与完整五量局部块试算，明确两者均未通过内部收敛。"""
from pathlib import Path
import os, json, hashlib

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
target = VAL / 'block_linear_research_summary_v7.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


out = read(VAL / 'block_linear_research_summary_v6.json')
for name, code in [('rowmomentumreal01', 0), ('coupledpatchreal01', 0), ('coupledpatchunit03', 0)]:
    process = read(VAL / (name+'_process.json'))
    assert process['status'] == 'terminal' and process['exitCode'] == code
    assert all(process['sourceStableAfter'].values())

row_file = VAL / 'rowmomentumreal01/row_momentum_real_summary.json'
row = read(row_file)
row_audit_file = VAL / 'rowmomentumreal01_independent_row_momentum.json'
row_audit = read(row_audit_file)
assert row_audit['status'].startswith('passed-') and row['numericalRunComplete']
assert not row['originalUnitTargetMet'] and not row['rowScaledTargetMet'] and not row['originalAResidualImproved']
assert row['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
row_inner = row['innerReport']
out['studies'].append(dict(label='原A三动量 / 等价行尺度', caseID='rowmomentumreal01',
    linearSystem='equivalent-R-A-with-P-R-inverse-right-preconditioning', arnoldiSteps=32,
    CApplications=0, ILUApplications=0, innerRelativeCResidual=None,
    basisOrthogonalityError=row_audit['basisOrthogonalityError'],
    modelRelationRelativeDifference=row_audit['rowSystemCombinationRelationRelativeDifference'],
    orthogonalityThreshold=1e-8, modelRelationThreshold=1e-8, candidateUsable=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, evidenceFile=str(row_file), evidenceSHA256=sha(row_file),
    displayState='原单位残差未改善',
    note='32步：尺度残差0.999779、原A残差1.00168，均未达1e-3；原A较未缩放0.999884变差。33次尺度A作用及1次原单位诊断另列，不是完整C或原NS收敛。'))
out['rowMomentumScaleDiagnostics'] = dict(auditFile=str(row_audit_file), auditSHA256=sha(row_audit_file),
    originalAResidual=row_audit['originalAResidual'], rowScaledResidual=row_audit['rowScaledResidual'],
    originalUnitTargetMet=False, rowScaledTargetMet=False, originalAResidualImproved=False,
    all32PhysicalArnoldiColumnsIndependentlyVerified=True, originalTargetNotReplacedByScaledFlag=True,
    factorReusedByByteIdentity=True, newFactorConstructions=0, newSchurAssemblies=0,
    solverAApplications=33, originalUnitDiagnosticAApplications=1, fixedLocalMapApplications=32,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0)

coupled_file = VAL / 'coupledpatchreal01/coupled_patch_real_summary.json'
coupled = read(coupled_file)
coupled_audit_file = VAL / 'coupledpatchreal01_independent_coupled_patches.json'
audit = read(coupled_audit_file)
assert audit['status'].startswith('passed-') and coupled['innerRunComplete'] and not coupled['candidateRejected']
assert not coupled['internalCTargetMet'] and not audit['internalCTargetMet']
assert coupled['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
inner = coupled['innerReport']
cycle = inner['krylov']['cycles'][0]
out['studies'].append(dict(label='完整C五变量 / 强连接局部块', caseID='coupledpatchreal01',
    linearSystem='C=-Jgraph-with-all-five-conservation-fields', arnoldiSteps=32,
    CApplications=33, ILUApplications=0, innerRelativeCResidual=audit['finalFullCResidual'],
    basisOrthogonalityError=cycle['basisOrthogonalityError'],
    modelRelationRelativeDifference=cycle['arnoldiRelationRelativeDifference'],
    orthogonalityThreshold=1e-8, modelRelationThreshold=1e-8, candidateUsable=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, evidenceFile=str(coupled_file), evidenceSHA256=sha(coupled_file),
    displayState='完整C内部未收敛',
    note='全部217272CV、29523个连通局部组，每组最多8CV及全部五量。32步原C残差0.986936，未达1e-3；全部局部逆/最终C已独立核对，未保存真实完整V/Z，不宣称逐列独立审计或原NS通过。'))
out['coupledPatchDiagnostics'] = dict(auditFile=str(coupled_audit_file), auditSHA256=sha(coupled_audit_file),
    cellCount=217272, scalarUnknowns=1086360, patchCount=29523, maximumPatchCells=8,
    allOriginalCVsExactlyOnce=True, allFiveFieldsCoupled=True, physicalGridChanged=False,
    allOriginalCConnectionsRetained=True, originalCVDiagonalFallbackPatches=0,
    maximumCertifiedOriginalInputInverseDefect=audit['maximumCertifiedOriginalInputInverseDefect'],
    maximumOrdinaryBLASInverseDefect=audit['maximumOrdinaryBLASInverseDefect'], decimal80CheckedPatches=audit['decimal80CheckedPatches'],
    originalCResidual=audit['finalFullCResidual'], originalCInternalTarget=1e-3, internalCTargetMet=False,
    solverCApplications=33, diagnosticCApplications=2, fixedLocalMapApplications=33,
    nativeArnoldiGatesPassed=True, realFullArnoldiBasisSaved=False, independentAllRealArnoldiColumnAuditClaimed=False,
    completeFiveFieldOrOriginalNSConvergenceProven=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0)
out['status'] = 'source-bound-row-scale-and-bounded-all-five-field-patches-with-unresolved-original-residuals'
out['diagnosticBuilderProcessID'] = os.getpid()
out['builderSHA256'] = sha(__file__)
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], researchRows=len(out['studies']), bytes=target.stat().st_size)))
