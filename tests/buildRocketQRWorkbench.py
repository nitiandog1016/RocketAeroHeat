"""将实际通过原NS的QR方向及其来源接回统一工作台，逐值核对使用方向。"""
from pathlib import Path
import json, shutil, os, sys
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, sha
OUT=ROOT/'outputs/v85_validation'


def alias():
    study=OUT/'qr_real01/augmented_qr_study.json';r=json.loads(study.read_text('utf8'))
    assert r['main'] and r['main']['acceptedNonlinearIterations']==1
    source=OUT/'qr_real01/krylov01';target=ROOT/'outputs/v83_validation/krylov26'
    assert not target.exists();target.mkdir()
    for name in ['wall_velocity_krylov_evidence.mat','wall_velocity_krylov_summary.json']:
        os.link(source/name,target/name);assert sha(source/name)==sha(target/name)
    result=dict(status='byte-identical-display-alias-of-actual-QR-accepted-state',sourceFolder=str(source),
                displayCase='krylov26',fileSHA256=sha(target/'wall_velocity_krylov_evidence.mat'),
                newFlowAcceptances=0,newProductionSolvedAreaM2=0)
    (OUT/'qr_view_source_binding.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result))


def bind(study_file,candidate_file,state_file):
    r=json.loads(study_file.read_text('utf8'))
    c=loadmat(candidate_file,simplify_cells=True,variable_names=['value','detail','preparationFile','preparationSHA256','sourceBinding'])
    s=loadmat(state_file,simplify_cells=True,variable_names=['result','report']);q=s['result']
    p=json.loads(Path(c['preparationFile']).read_text(encoding='utf-8-sig'))
    assert sha(c['preparationFile'])==c['preparationSHA256']==r['preparationSHA256']
    assert c['sourceBinding']==r['sourceBinding']==p['sourceBinding']
    assert sha(p['sourceBinding']['sourceFile'])==p['sourceBinding']['sourceFileSHA256']
    lines=q['linearLog'];line=lines if isinstance(lines,dict) else lines[0]
    np.testing.assert_array_equal(c['value'],line['steadyDirectionEvidence']['directionU'])
    probe=loadmat(q['options']['steadyInitialGuess']['file'],simplify_cells=True,variable_names=['sourceBinding','directionU'])
    np.testing.assert_array_equal(c['value'],probe['directionU'])
    assert probe['sourceBinding']['sha256']==c['sourceBinding']['sourceFileSHA256']
    assert r['modelPointOneMet'] and r['main']['residualEvaluations']==s['report']['residualEvaluations']
    detail=c['detail'];assert detail['method']=='complete-parent-C-fixed-primitive-penalty-augmented-patch-QR-LSQR'
    return dict(kind='augmented-patch-QR-LSQR-v1',methodLabel=f"完整增广像分块QR右预处理；零初猜{detail['nativeIterations']}步LSQR",
                studyFile=str(study_file),studyFileSHA256=sha(study_file),candidateFile=str(candidate_file),
                candidateFileSHA256=sha(candidate_file),baseCase='krylov25',usedTotalDirectionExactlyVerified=True,
                hardUpdateLimitClaimed=False,firstSegmentIterations=detail['nativeIterations'],restartedSegmentIterations=0,
                internalCRelativeResidual=detail['relativeCResidual'],internalCTarget=detail['options']['relativeTolerance'],
                internalCTargetMet=bool(detail['internalCTargetMet']),completeModelRelativeDefect=r['completeSelectedModelRelativeDefect'],
                originalNSDirectionDefect=line['fullSystemRelativeDefect'],originalNSDirectionTarget=q['options']['linearRelativeTolerance'])


def build():
    folder=OUT/'workbench01';source=OUT/'workbench_base';assert not folder.exists()
    prep=bind(OUT/'qr_real01/augmented_qr_study.json',OUT/'qr_real01/augmented_qr_direction.mat',
              ROOT/'outputs/v83_validation/krylov26/wall_velocity_krylov_evidence.mat')
    shutil.copytree(source,folder)
    d=json.loads((folder/'wall-velocity-study-data.json').read_text('utf8'))
    d['cases'][2]['label']='v85 分块QR实际接受态（未收敛）'
    d['cases'][2]['solverDiagnostics']['directionPreparation']=prep
    raw=json.dumps(d,ensure_ascii=False,separators=(',',':'))
    (folder/'wall-velocity-study-data.json').write_text(raw,encoding='utf8')
    (folder/'wall-velocity-study-data.js').write_bytes(('window.ROCKET_WALL_VELOCITY_STUDY='+raw+';\n').encode('utf8'))
    identical=[]
    for file in sorted((source/'cases').glob('*/results.json')):
        rel=file.relative_to(source);assert file.read_bytes()==(folder/rel).read_bytes();identical.append(str(rel))
    receipt=dict(status='bound-actual-QR-direction-and-original-NS-state-in-five-step-workbench',
                 directionPreparation=prep,originalFacetDatasetsByteIdentical=identical,
                 physicalSolution=False,newProductionSolvedAreaM2=0,browserVisualOrWebGLRuntimeVerified=False)
    (folder/'qr_view_binding.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(receipt,ensure_ascii=False))


if __name__=='__main__':
    {'alias':alias,'build':build}[sys.argv[1]]()
