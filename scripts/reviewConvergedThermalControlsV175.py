"""只读v105/v108准确小回执和源，核对旧收敛冷壁控制适用范围；不重算或读体场。"""
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v175_validation/control_review02'
PLAN=ROOT/'docs/validation/v175/PLAN.md'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def config(text):
    result={}
    for line in text.splitlines():
        line=line.split('%')[0].strip()
        if not line or '=' not in line:continue
        key,value=[v.strip() for v in line.split('=',1)];assert key not in result
        if re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?',value):value=float(value)
        result[key]=value
    return result


def main():
    start=time.monotonic();assert not OUT.exists();OUT.mkdir(parents=True)
    shutil.copy2(str(Path(__file__)),str(OUT/'executed_review.py'))
    locks={};report=dict(passed=False,newCFDRuns=0,newBuilds=0,largeFieldReads=0,
        oldRegressionReruns=0,similaritySolverReruns=0,newPhysicalQualification=False)
    def load(path,expected=None,as_json=True):
        path=Path(path);assert path.stat().st_size<2*1024**2,'只允许明确小证据，不读大场'
        value=sha(path)
        if expected is not None:assert value==expected,'SHA不匹配: '+str(path)
        locks[str(path.resolve())]=value
        text=path.read_text(encoding='utf-8-sig')
        return json.loads(text) if as_json else text
    try:
        load(Path(__file__),as_json=False);load(PLAN,as_json=False)
        manifests={v:load(ROOT/('docs/validation/%s/delivery_manifest.json'%v)) for v in ('v105','v108')}
        evidence={v:{row['path']:row for row in manifests[v]['evidence']} for v in manifests}
        def archived(v,relative,as_json=True):
            row=evidence[v][relative]
            return load(ROOT/relative,row['sha256'],as_json)
        def source(v,prefix):
            matches=[(digest,name) for digest,name in manifests[v]['sourceVersions'].items()
                     if Path(name).name.startswith(prefix)]
            assert len(matches)==1,prefix
            digest,name=matches[0];text=load(ROOT/name,digest,False)
            return dict(path=name,SHA256=digest),text
        coarse_cfg=archived('v105','docs/validation/v105/cases/plate03/case.cfg',False)
        fine_cfg=archived('v108','docs/validation/v108/engines/warm/case.cfg',False)
        a=config(coarse_cfg);b=config(fine_cfg)
        coarse=archived('v105','docs/validation/v105/cases/plate03/case_manifest.json')
        fine=archived('v108','docs/validation/v108/engines/warm/case_manifest.json')
        cr=archived('v105','docs/validation/v105/cases/plate03/engine_process.json')
        fr=archived('v108','docs/validation/v108/engines/warm/engine_process.json')
        ca=archived('v105','docs/validation/v105/plate03_audit_final.json')
        fa=archived('v108','docs/validation/v108/checks/plate_warm01/physical_audit.json')
        warm=archived('v108','docs/validation/v108/checks/plate_warm01/warm_start_manifest.json')
        comparison=archived('v108','docs/validation/v108/checks/plate_warm01/mesh_comparison.json')
        init_audit=archived('v108','docs/validation/v108/checks/plate_warm01/initial_audit.json')
        builder,builder_text=source('v105','buildSU2FlatPlateCase.m.11359cb40830a638')
        initializer,init_text=source('v108','prepareSU2FlatPlateWarmStart.py.4865173ac3120306')
        old_audit,audit_text=source('v105','auditEmbeddedSU2FlatPlate.py.9dfc4d744b648be4')
        new_audit,new_audit_text=source('v108','auditEmbeddedSU2FlatPlate.py.32c0d70089bf4d2a')
        warm_check,warm_text=source('v108','auditSU2WarmConvergence.py.8f30d49d1137354b')
        assert ca['auditSHA256']==old_audit['SHA256'] and fa['auditSHA256']==new_audit['SHA256']
        assert initializer['SHA256']==warm['inputSHA256'][str(ROOT/'scripts/prepareSU2FlatPlateWarmStart.py')]
        upstream={}
        for name in ('CNSSolver.cpp','CFVMFlowSolverBase.inl','flatplate_supersonic.cfg'):
            relative='docs/validation/v105/upstream/'+name
            upstream[name]=archived('v105',relative,False)
        assert 'Apply a weak boundary condition for the energy equation' in upstream['CNSSolver.cpp']
        assert 'thermal_conductivity * (There - Twall_local) / dist_ij * RefHeatFlux' in upstream['CFVMFlowSolverBase.inl']
        assert '(rest[:, 6] / rho - (u*u+v*v)/2) * (gamma-1)/R' in audit_text
        assert 'strictlyConservativeRemap=False' in init_text and '物理/离散配置须一致' in init_text
        # 当前文件只用作同名历史本地输入的身份核对；不调用当前准备/审核实现。
        for folder,receipt,cfg,metadata in [('outputs/v105_validation/plate03',cr,coarse_cfg,coarse),
                                           ('outputs/v108_validation/plate_warm01',fr,fine_cfg,fine)]:
            folder=ROOT/folder
            assert load(folder/'engine_process.json')==receipt
            assert load(folder/'case.cfg',as_json=False)==cfg
            assert load(folder/'case_manifest.json')==metadata
            assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut']
            assert all(receipt['inputsStable'].values())
            assert receipt['inputSHA256'][str(folder/'case.cfg')]==sha(folder/'case.cfg')
        ignored={'ITER','OUTPUT_WRT_FREQ','WRT_RESTART_OVERWRITE','WRT_SURFACE_OVERWRITE',
                 'RESTART_SOL','SOLUTION_FILENAME','READ_BINARY_RESTART'}
        assert {k:v for k,v in a.items() if k not in ignored}=={k:v for k,v in b.items() if k not in ignored}
        assert a['RESTART_SOL']=='NO' and b['RESTART_SOL']=='YES' and b['SOLUTION_FILENAME']=='initial_flow'
        assert a['FLUID_MODEL']=='IDEAL_GAS' and a['KIND_TURB_MODEL']=='NONE'
        assert ca['residualConverged'] and fa['residualConverged']
        for review in (ca,fa):
            assert review['lastHistoryFields']['rmsRho']<=-8 and review['lastHistoryFields']['rmsRhoE']<=-8
            assert not review['productionQualified'] and not review['fullNavierStokesVerified']
        assert not comparison['gridIndependenceQualified'] and not comparison['experimentCompared']
        assert comparison['pressureScreenPassed'] and comparison['heatFluxScreenPassed']
        gamma=a['GAMMA_VALUE'];gas=a['GAS_CONSTANT'];cv=gas/(gamma-1);cp=gamma*cv
        temperature=a['FREESTREAM_TEMPERATURE'];mach=a['MACH_NUMBER'];total=temperature*(1+(gamma-1)*mach**2/2)
        recovery=temperature+math.sqrt(a['PRANDTL_LAM'])*(total-temperature)
        assert math.isclose(total,1260.,rel_tol=1e-15) and 300<recovery<total
        actual_manifest=load(ROOT/'outputs/v174_validation/real01/case_manifest.json')
        actual_cfg=config(load(ROOT/'outputs/v174_validation/real01/case.cfg',as_json=False))
        actual_review=load(ROOT/'outputs/v174_validation/independent_real01/report.json')
        assert actual_review['passed'] and actual_review['termination']['modelDomainExceeded']
        changes={key:dict(control=a.get(key),v174=actual_cfg.get(key)) for key in sorted(set(a)|set(actual_cfg))
                 if a.get(key)!=actual_cfg.get(key)}
        report.update(passed=True,exactHistoricalSourceReview=True,sourceBindings=[builder,initializer,old_audit,new_audit,warm_check],
            nativeRuntime=dict(version=cr['executableVersion'],archiveSHA256=cr['archiveSHA256'],
                executableSHA256=cr['inputSHA256'][str(ROOT/'tools/runtime/su2/8.6.0/bin/SU2_CFD.exe')],
                currentExecutableRehashed=False,meaning='original signed-by-hash execution receipt; no runtime rebuild/retest'),
            controls=[dict(name='v105/plate03',configuration=a,meshMetadata=coarse,
                nativeSeconds=cr['elapsedSeconds'],lastHistory=ca['lastHistoryFields'],heatReference=ca['reference'],
                referenceExtended=ca['referenceExtended'],comparisonInterval=ca['comparisonInterval'],
                boundaryLayerDifferenceRMS=ca['boundaryLayerRelativeDifferenceRMS'],leadingEdge=ca['leadingEdge'],
                wallNodeTemperatureDeviationMaxK=ca['fluidNodeMinusPrescribedTemperatureMax'],
                interiorTemperatureDeviationMaxK=ca['fluidNodeMinusPrescribedTemperatureInteriorMax']),
              dict(name='v108/plate_warm01',configuration=b,meshMetadata=fine,nativeSeconds=fr['elapsedSeconds'],
                lastHistory=fa['lastHistoryFields'],heatReference=fa['reference'],leadingEdge=fa['leadingEdge'],
                boundaryLayerDifferenceRMS=fa['boundaryLayerRelativeDifferenceRMS'],
                wallNodeTemperatureDeviationMaxK=fa['fluidNodeMinusPrescribedTemperatureMax'],
                interiorTemperatureDeviationMaxK=fa['fluidNodeMinusPrescribedTemperatureInteriorMax'])],
            energyAndWall=dict(reference='native calorically perfect e=cv*T, h=cp*T, no OTH offset',cvJkgK=cv,cpJkgK=cp,
                freestreamT=temperature,prescribedTw=300.,stagnationTemperatureK=total,
                laminarFlatPlateRecoveryEstimateK=recovery,referenceHeatSign='positive into cold wall',
                controlIsHeatedAboveFreestreamStagnation=False,wallTemperatureIsWeakFluxInput=True,
                outputTemperatureIsSolvedFluidNode=True,wallFormula='q_into_wall=k_wall*(T_normal_neighbor-Tw)/normal_distance'),
            warmStart=dict(method=warm['method'],sourceNodes=warm['sourceNodes'],targetNodes=warm['targetNodes'],
                integralRelativeChanges=warm['integralRelativeChanges'],strictlyConservativeRemap=warm['strictlyConservativeRemap'],
                minTemperatureK=warm['minTemperature'],maxTemperatureK=warm['maxTemperature'],
                samePhysicsAndNumericsBeforeRestart=True,requiresConvergedSameConditionCoarseSource=True,
                noTemperatureClipping=True,noFieldExtrapolation=True,newWallBoundaryFormula=False),
            numericalReference=dict(qSqrtX=ca['reference']['qSqrtX'],qLaw='q_ref(s)=12264.843223508966/sqrt(s) W/m^2, s in m',
                comparisonIntervalM=ca['comparisonInterval'],twoGrid=comparison['changes'],
                experiment=False,independentFullNavierStokesTruth=False,gridIndependence=False,
                leadingEdgeQualified=False,completeSurfaceQualified=False,newRocketQualifiedAreaM2=0),
            v174Comparison=dict(configurationDifferences=changes,
                givenTwRangeK=actual_manifest['publicBenchmark']['conditions']['wallTemperatureRangeK'],
                stagnationTemperatureK=actual_manifest['publicBenchmark']['conditions']['totalTemperature_K'],
                actualDomainFailurePreserved=True,sameCaloricEnergyZeroConvention=True,
                conservativeStateDirectlyReusable=False,geometryAndBoundaryLayerSourceMissing=True,
                oldControlsDoNotAuthorizeLargerV174DevelopmentDomain=True),
            decision=dict(portableCompleteHeatedNACAStartupEstablished=False,
                establishedCapability='same-geometry and same-condition converged coarse-to-fine cold-wall control',
                actionableNewNACACandidateFromTheseControls=False,
                reason='different thermal direction and geometry; same weak boundary; no converged matching NACA seed',
                fluxOrCFLScanRecommended=False,newCFDAuthorized=False),
            exclusions=dict(v129='different SST-2003m/M6Tw025 DNS comparison; not these converged laminar controls',
                largeSavedFieldsNotRevalidated=True,meshGeometryNotRecomputed=True,oldReportsNotRewritten=True))
        assert time.monotonic()-start<40
    finally:
        report['inputSHA256']=locks;report['inputsStable']={p:sha(p)==digest for p,digest in locks.items()}
        report['passed']=report['passed'] and all(report['inputsStable'].values())
        report['elapsedSeconds']=time.monotonic()-start
        (OUT/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
        print(json.dumps(dict(passed=report['passed'],inputs=len(locks),elapsedSeconds=report['elapsedSeconds']),ensure_ascii=False))


if __name__=='__main__':main()
