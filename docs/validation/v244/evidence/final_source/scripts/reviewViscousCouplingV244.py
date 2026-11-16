"""模块职责：选择保留真实下游的成熟内部耦合，绑定下一整体控制的准确来源。"""
from pathlib import Path
import sys
from runSU2RocketPhaseV222 import ROOT,read,save,sha
import prepareViscousRegionV244 as run
BASE=run.BASE;DOC=run.DOC;CLOSED=True
OLD=ROOT/'tools/runtime/su2-thermal-build/source'
SPECS=[
 ('SU2_CFD/include/solvers/CScalarSolver.inl','void CScalarSolver<VariableType>::BC_Riemann',20),
 ('SU2_CFD/src/solvers/CSolver.cpp','void CSolver::GetCommCountAndType',49),
 ('Common/src/geometry/CPhysicalGeometry.cpp','std::unique_ptr<CADTElemClass> CPhysicalGeometry::ComputeViscousWallADT',89),
 ('Common/src/adt/CADTElemClass.cpp','if (globalTree)',84)]
MG=[('void CMultiGridIntegration::MultiGrid_Cycle',118),
    ('void CMultiGridIntegration::SetForcing_Term',56),
    ('void CMultiGridIntegration::SetRestricted_Solution',58)]


def review():
    """接口入射数据缺失时保留完整物理域；并行分区不冒称收敛修复。"""
    out=BASE/'coupling01';out.mkdir(exist_ok=False)
    audited=read(BASE/'audit02/report.json');assert audited['passed']
    known=read(ROOT/'outputs/v167_validation/build01/source_tree.json');rows=[]
    specs=[]
    for rel,anchor,count in SPECS:
        expected=known[str(OLD/rel)]
        candidates=[OLD/rel,ROOT/'outputs/v167_validation/build01/executed_sources'/rel]
        p=next((p for p in candidates if p.is_file() and sha(p)==expected),None)
        assert p,('EXACT_SOURCE_MISSING',rel);specs.append((p,anchor,count,expected))
    rel='SU2_CFD/src/integration/CMultiGridIntegration.cpp'
    bound=next(r for r in read(run.SOURCES/'bindings.json') if r['relative']==rel)
    specs += [(run.SOURCES/rel,anchor,count,bound['sha256']) for anchor,count in MG]
    for k,(p,anchor,count,expected) in enumerate(specs):
        assert sha(p)==expected;lines=p.read_text(encoding='utf8').splitlines()
        positions=[i for i,s in enumerate(lines) if anchor in s];assert len(positions)==1
        first=positions[0];last=min(first+count,len(lines));dest=out/('%02d_%s.txt'%(k,p.name))
        dest.write_text('\n'.join('%6d %s'%(i+1,lines[i]) for i in range(first,last))+'\n',encoding='utf8')
        rows.append(dict(path=str(p),sha256=expected,firstLine=first+1,lastLine=last,
            excerpt=dest.name,excerptSHA256=sha(dest)))
    mesh=ROOT/'outputs/v167_validation/real01/mesh.su2'
    expected=read(ROOT/'docs/validation/v243/evidence/final_report.json')['assets'][mesh.relative_to(ROOT).as_posix()]
    assert sha(mesh)==expected
    contract=dict(schema='RAH.ConnectedViscousProvider.v1',geometry=str(mesh),meshSHA256=expected,
        physicalZones=1,originalPhysicalWalls=19752,originalCells=705490,originalNodes=335726,
        upstreamCells=160304,downstreamCells=545186,internalInterfaceFaces=1984,
        decompositionIndex=str(run.INDEX),indexSHA256=sha(run.INDEX),
        partitionIsObservationOnly=True,artificialPhysicalMarkers=[],
        downstreamFlowIsSolvedInternally=True,externalCFDImportRequired=False,
        freestreamAndWallInputsFrom=str(ROOT/'outputs/v222_validation/rocket01/case.cfg'),
        stateInitialization='compatible uniform freestream in any separately admitted control',
        old300StateUsedAsBoundary=False,old300RestartAdmitted=False,
        originalNearWallAssetsChanged=False,sourceContract='phase-Air / physical energy offset 300000 / SA-neg WITHFT2',
        modelCombination='same native single-zone interior flow-SA operators; MZ17 not used',
        fullDomainPeakReferenceBytes=2153095168,referenceIsNotNewPeakGuarantee=True,
        currentRunAdmitted=False,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'connected_provider.json',contract)
    decision=dict(passed=True,sourceExcerpts=rows,physicalCutRejected=True,
        missingDownstreamFaceCount=audited['incomingUpstreamCounts']['1']+audited['incomingUpstreamCounts']['4']+audited['incomingUpstreamCounts']['5'],
        characteristicStaticPressureExistsInFlowCode=True,phase15StaticPressureConfigurationAllowed=False,
        SAStaticPressureDispatch='BC_Outlet with own scalar copy; no diffusive boundary flux',
        viscousRiemannSecondarySupported=True,
        adoptedPhysicalClosure='retain downstream in the same native physical zone; cut becomes observation index',
        nativeMPICommunicationExists=True,nativeGlobalWallDistanceADTExists=True,
        mpiInstalledOrBuilt=False,mpiSelectedAsConvergenceFix=False,
        nextConvergenceCandidate='native FAS agglomeration with original mesh as finest level',
        fasRestrictionAndForcingSourceBound=True,fasPhaseAirActualCompatibilityTested=False,
        fasAdoptedForProduction=False,agglomerationDoesNotReplaceOriginalFineMesh=True,
        newBuilds=0,nativeQueries=0,flowUpdates=0,newQualifiedAreaM2=0)
    save(out/'report.json',decision)
    print(dict(passed=True,excerpts=len(rows),productionAdmitted=False,flowUpdates=0))


def main():
    if CLOSED:raise SystemExit('v244 source decision closed; use NEXT_DECISION')
    run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC;size=run.guard.size_bytes
    run.guard.size_bytes=lambda p:size(BASE) if Path(p)==BASE or Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
        ROOT/'outputs/v167_validation/build01/source_tree.json',run.SOURCES/'bindings.json',
        run.SOURCES/'SU2_CFD/src/integration/CMultiGridIntegration.cpp',BASE/'audit02/report.json',
        ROOT/'docs/validation/v243/evidence/final_report.json',run.INDEX,
        ROOT/'outputs/v167_validation/real01/mesh.su2',ROOT/'outputs/v222_validation/rocket01/case.cfg']
    for rel,_,_ in SPECS:
        locks += [p for p in (OLD/rel,ROOT/'outputs/v167_validation/build01/executed_sources'/rel) if p.is_file()]
    return run.guard.supervise('coupling',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,locks,BASE,'offline')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v244 source decision closed; use NEXT_DECISION')
    if sys.argv[1]=='task':review()
    else:raise SystemExit(0 if main() else 1)
