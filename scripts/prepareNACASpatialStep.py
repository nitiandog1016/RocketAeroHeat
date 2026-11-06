"""模块用途：同一已审核层流场接独立空间Tw限步包，保留原自适应参数与冷重启身份。"""
import argparse
import copy
import os
from pathlib import Path
import time
import numpy as np

import prepareNACAAdaptiveContinuation as original
from prepareSU2SpatialStepPatch import FLAVOR
from prepareSU2NACALaminar import capacity
from su2IdealAirInterface import CONSERVATIVE_STEP
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

KIND='naca-tn3344-spatial-step-cold-restart-v1'
INITIALIZATION='same-mesh-spatial-step-laminar-state-v1'
PLAN=ROOT/'docs/validation/v173/PLAN.md'
PROOF=original.PROOF
PROOF_SHA256=original.PROOF_SHA256
SPATIAL=original.SPATIAL
ADAPTIVE=copy.deepcopy(original.ADAPTIVE)
INITIAL_COLUMNS=original.INITIAL_COLUMNS
ACTUAL_PLAN_MARKER='SPATIAL_STEP_ACTUAL_CANDIDATE_REGISTERED_V173'
VARIANTS=('small-step','validation-step')
QUERY_CONTEXT=dict(kind='ideal-air-domain-query-context-v1',domainErrorTypePreserved=True,
    modes=['rhoe','PT','Prho','rhoT','hs','Ps'],legalStateUnchanged=True,
    scopedPaths=['CNSVariable.SetPrimVar','CEulerSolver.BC_Riemann'],unscopedExplicit=True,
    boundaryNumericsChanged=False,reconstructionPairRecoveryChanged=False)
NUMERICAL_CHANGE=dict(kind='laminar-node-conservative-backtracking-with-spatial-wall-v1',
    conservativeStep=copy.deepcopy(CONSERVATIVE_STEP),queryContext=copy.deepcopy(QUERY_CONTEXT),
    turbulence='NONE',turbulentK=0,eddyViscosity=0,
    runtimeParentFlavor=original.FLAVOR,originalCFLFeedbackChanged=False,
    effectiveBacktrackingFedToCFL=False,riemannBoundaryRecovery=False,
    wholeQueryPathAdmissibilityGuaranteed=False,physicalQualified=False)


def specification(variant):
    """模块一：只允许既有来源上的两个具名任务；实际候选需单独登记。"""
    assert variant in VARIANTS
    actual=variant=='validation-step'
    spec=original.specification('validation-adaptive' if actual else 'small-adaptive')
    spec.update(variant=variant,requiredRuntimeFlavor=FLAVOR,
        resourceProfile='v173-naca-spatial-step-validation' if actual else 'v173-naca-spatial-step-small')
    return spec


def source_state(variant):
    """模块二：复用准确v171来源证明，不从v172未知末态或未保存历史初始化。"""
    spec=specification(variant)
    result=original.source_state('validation-adaptive' if spec['actualValidationMesh'] else 'small-adaptive')
    result['spec']=spec
    return result


def source_files():
    names=('prepareNACASpatialStep.py','su2NACASpatialStepInterface.py','prepareSU2SpatialStepPatch.py',
        'prepareSU2ConservativeStepPatch.py')
    # 旧来源准备器及其原计划仍是独立来源证明的一部分。
    return original.source_files()+[ROOT/'scripts'/name for name in names]+[PLAN,
        ROOT/'third_party/ideal_air/ConservativeStep.inc',ROOT/'third_party/ideal_air/ConservativeStepTrace.hpp',
        ROOT/'third_party/ideal_air/PropertyQueryContext.hpp']


def adaptive_config(text,spec):
    return original.adaptive_config(text,spec)


def prepare(folder,variant):
    """模块三：独占目录、同网格只读硬链接与17位守恒初场，不改来源文件。"""
    from runEmbeddedSU2 import atomic_write_json
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs');folder.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();marker=folder/'step_preparation.json'
    atomic_write_json(marker,dict(status='preparing',exitCode=None))
    (folder/'executed_prepare.py').write_bytes(Path(__file__).read_bytes())
    (folder/'plan_at_prepare.md').write_bytes(PLAN.read_bytes())
    try:
        checked=source_state(variant);spec=checked['spec'];source=checked['source'];old=checked['manifest']
        if spec['actualValidationMesh']:
            assert ACTUAL_PLAN_MARKER in PLAN.read_text(encoding='utf8'),'实际候选尚未事前登记'
        locks=dict(checked['locks']);locks.update({str(p):sha(p) for p in source_files()})
        for name in ('mesh.su2','mesh_arrays.npz','geometry.json','wall_temperature.csv'):
            os.link(source/name,folder/name)
        np.savetxt(str(folder/'initial_flow.csv'),checked['initial'],delimiter=',',fmt='%.17g',comments='',
            header=','.join('"'+name+'"' for name in INITIAL_COLUMNS))
        initial=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1)
        assert np.array_equal(initial,checked['initial'])
        (folder/'case.cfg').write_text(adaptive_config(checked['text'],spec),encoding='ascii')
        initialization=dict(checked['identity'],kind=INITIALIZATION,localIterationStartsAt=0,
            initialStateSHA256=sha(folder/'initial_flow.csv'),sourceSHA256=locks,
            samePhysicalModel=True,sameSpatialDiscretization=True,samePhysicsAndNumerics=False,
            numericalChange=copy.deepcopy(NUMERICAL_CHANGE),uninterruptedTrajectoryEquivalent=False,
            adaptiveCFLReset=True,solverHistoryReset=True,boundaryLayerSolved=False,physicsQualified=False,newQualifiedArea=0)
        save(folder/'step_initialization.json',initialization)
        wallids=np.unique([i for face in checked['markers']['wall'] for i in face[1:]])
        contract=copy.deepcopy(checked['contract'])
        contract.update(schemaVersion=8,initialization=INITIALIZATION,
            initializationSHA256=sha(folder/'step_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
            checkpointContract=capacity(len(initial),len(wallids),len(checked['markers']['wall']),
                spec['updates'],spec['every'],spec['maximumOutputBytes']))
        save(folder/'ideal_air_contract.json',contract)
        manifest=dict(caseKind=KIND,specification=spec,requiredRuntimeFlavor=FLAVOR,
            idealAir=contract,idealAirContractSHA256=sha(folder/'ideal_air_contract.json'),
            initialization=initialization,advancement=copy.deepcopy(ADAPTIVE),
            numericalChange=copy.deepcopy(NUMERICAL_CHANGE),
            thermalBoundary=copy.deepcopy(old['thermalBoundary']),thermalProfileSHA256=old['thermalProfileSHA256'],
            publicBenchmark=copy.deepcopy(old['publicBenchmark']),nodes=old['nodes'],elements=old['elements'],
            expectedVelocity=old['expectedVelocity'],physicalQualified=False,planSHA256=sha(PLAN),
            historicalProofSHA256=PROOF_SHA256)
        save(folder/'case_manifest.json',manifest)
        request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
            maximumSeconds=spec['maximumSeconds'],threads=spec['threads'],resourceProfile=spec['resourceProfile'],
            sourceFiles=[p.relative_to(ROOT).as_posix() for p in source_files()])
        save(folder/'request.json',request)
        assert all(sha(Path(p))==h for p,h in locks.items())
        atomic_write_json(marker,dict(status='terminal',exitCode=0,elapsedSeconds=time.monotonic()-started,
            initializationSHA256=contract['initializationSHA256'],sourceInputsStable=True))
        return manifest
    except Exception as error:
        atomic_write_json(marker,dict(status='terminal',exitCode=1,elapsedSeconds=time.monotonic()-started,
            error=type(error).__name__+': '+str(error)))
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True)
    parser.add_argument('--variant',choices=VARIANTS,required=True);args=parser.parse_args()
    prepare(args.folder,args.variant)
