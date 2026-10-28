"""独立审查双温域选择、冷却释放、活跃循环、完整增广目标与稠密最终子空间。"""
from pathlib import Path
from datetime import datetime,timezone
import os,sys,json,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.linalg import null_space
from screenRocketFullFLeastSquares import ROOT,sha,local
from screenRocketParentImageScale import setup
from screenRocketAugmentedPatchQR import metric_blocks,block_weight
from auditRocketTemperatureProjection import row_reference
from rocketTemperatureProjectorReference import project_coordinates,certify_projector_span


def cells(x):
    if isinstance(x,dict):return [x]
    return list(np.atleast_1d(x).ravel()) if np.size(x) else []


def active_arithmetic(rows,evidence,report):
    H=rows['H'];signs=np.atleast_1d(rows['temperatureBoundSigns']);free=evidence['free'];active=np.array([],dtype=int)
    def inspect(v):
        t=np.atleast_1d(H@v.ravel(order='F'));b=128*np.finfo(float).eps*np.atleast_1d(abs(H)@abs(v.ravel(order='F')))
        return t,b,np.flatnonzero(signs*t>b)
    tangent,bound,violating=inspect(free);initial_violating=(violating+1).tolist();counts=int(evidence['freeDetail']['nativeIterations'])
    cycle_checks=[]
    for cycle in cells(evidence['cycles']):
        expected=np.unique(np.r_[active,violating]);actual=np.atleast_1d(cycle['activeRowIndices']).astype(int)-1
        np.testing.assert_array_equal(actual,expected);active=actual
        tangent,bound,violating=inspect(cycle['value']);np.testing.assert_allclose(tangent,np.atleast_1d(cycle['temperatureDerivativeK']),rtol=1e-8,atol=1e-8)
        np.testing.assert_allclose(bound,np.atleast_1d(cycle['localDotRoundoffBoundK']),rtol=1e-12,atol=0)
        np.testing.assert_array_equal(violating+1,np.atleast_1d(cycle['remainingViolatingRowIndices']).astype(int))
        projector=cycle['projector'];exact_identity=None
        if 'originalToUniqueRow' in projector:
            representatives=np.atleast_1d(projector['exactRepresentativeRowIndices']).astype(int)-1
            mapping=np.atleast_1d(projector['originalToUniqueRow']).astype(int)-1
            selected=H[active];assert (selected-selected[representatives[mapping]]).nnz==0
            assert projector['constraintCount']==len(active) and projector['independentConstraintCount']==len(representatives)
            assert projector['allOriginalRowsExactlyRepresented'] and not projector['rankTruncationUsed']
            exact_identity=dict(originalRows=len(active),independentConditions=len(representatives),
                                allOriginalRowsExactlyReconstructed=True,rankTruncationUsed=False)
        counts+=int(cycle['detail']['nativeIterations']);cycle_checks.append(dict(activeRows=(active+1).tolist(),remainingViolatingRows=(violating+1).tolist(),exactDuplicateRepresentation=exact_identity))
    assert report['bandDirectionFeasible']==(len(violating)==0) and report['totalNativeIterations']==counts
    np.testing.assert_allclose(tangent,np.atleast_1d(report['temperatureDerivativeK']),rtol=1e-8,atol=1e-8)
    assert not report['oldActiveSetInherited'] and not report['constrainedOptimumClaimed'] and not report['physicalBoundaryConditionAdded']
    return dict(initialOutwardRows=initial_violating,cycles=cycle_checks,totalNativeIterations=counts,allCurrentBandRowsFeasible=len(violating)==0)


def unit(folder,read):
    s=read(folder/'temperature_active_unit.mat');checks=[]
    public=row_reference(s['bandRows'])
    assert set(np.atleast_1d(s['bandRows']['temperatureBoundSigns']).tolist())=={-1,1}
    for c in cells(s['cases']):
        row=c['rows'];r=c['detail'];active=np.atleast_1d(r['activeRowIndices']).astype(int)-1
        H=row['H'][active].toarray();b=c['source'];penalty=c['factor']['penaltyInOriginalResidualUnits'];rhs=b['base']['rhs'].ravel(order='F');N=len(rhs)
        A=np.row_stack([b['base']['C'],penalty*b['W']]);Z=null_space(H);dense=Z@np.linalg.lstsq(A@Z,np.r_[rhs,np.zeros(N)],rcond=None)[0]
        error=float(np.linalg.norm(c['value'].ravel(order='F')-dense)/max(1.,np.linalg.norm(dense)));assert error<1e-8
        detail=active_arithmetic(row,c['evidence'],r);assert 2 not in detail['cycles'][0]['activeRows']
        checks.append(dict(denseFinalActiveSubspaceDifference=error,activeArithmetic=detail))
    return dict(publicTemperatureRows=public,cases=checks,inequalityOptimumClaimed=False)


def real(folder,read,study_file=None):
    study=json.loads((study_file or folder/'active_continuation.json').read_text('utf8'));checks=[]
    for step in study['steps']:
        sf=Path(step['directionFile']).parent;s=read(step['directionFile']);context=read(sf/'direction_context.mat',['factor'])
        p=json.loads(Path(s['preparationFile']).read_text(encoding='utf8'));assert sha(s['preparationFile'])==s['preparationSHA256']==step['preparationSHA256']
        assert sha(p['sourceBinding']['sourceFile'])==p['sourceBinding']['sourceFileSHA256']==s['sourceBinding']['sourceFileSHA256']
        for item in p['artifacts'].values():assert sha(item['file'])==item['sha256']
        packed=read(p['artifacts']['packed']['file'],['packed'])['packed'];factor=context['factor'];bs=cells(factor['buckets'])
        for b in bs:b['inverse']=b['map']
        factor['buckets']=bs;A,plan,order=setup(packed,factor);del packed,factor,bs;gc.collect()
        n=len(order);metric=metric_blocks(s['metric'],order);rhs=s['rhs'][order].ravel();penalty=context['factor']['penaltyInOriginalResidualUnits']
        data=[('free',s['evidence']['free'],s['evidence']['freeCoordinates'],s['evidence']['freeDetail'],None)]
        for j,c in enumerate(cells(s['evidence']['cycles'])):data.append((f'active{j+1}',c['value'],c['coordinates'],c['detail'],c['projector']))
        directions=[]
        for label,value,coordinates,detail,projector in data:
            projected=coordinates if projector is None else project_coordinates(projector,coordinates)
            y=projected[order].ravel()
            lift=local(plan,y);v=value[order].ravel();lift_error=float(np.linalg.norm(lift-v)/max(1.,np.linalg.norm(v)));assert lift_error<1e-12
            image=A@v;w=block_weight(metric,v);relative=float(np.linalg.norm(rhs-image)/np.linalg.norm(rhs));objective=relative**2+(penalty*np.linalg.norm(w)/np.linalg.norm(rhs))**2
            np.testing.assert_allclose(relative,detail['relativeCResidual'],rtol=1e-10);np.testing.assert_allclose(objective,detail['normalizedAugmentedObjective'],rtol=1e-10)
            directions.append(dict(label=label,coordinateLiftDifference=lift_error,originalCResidual=relative,normalizedAugmentedObjective=objective))
        def transpose_physical(raw):
            algebraic=raw.reshape(n,5,order='F')[order].ravel();mapped=local(plan,algebraic,True)
            physical=np.empty((n,5));physical[order]=mapped.reshape(n,5);return physical.ravel(order='F')
        span_checks=[]
        for cycle in cells(s['evidence']['cycles']):
            active=np.atleast_1d(cycle['activeRowIndices']).astype(int)-1
            span_checks.append(certify_projector_span(cycle['projector'],s['rows']['H'][active],transpose_physical))
        arithmetic=active_arithmetic(s['rows'],s['evidence'],s['directionReport']);public=row_reference(s['rows']) if len(np.atleast_1d(s['rows']['faceIDs'])) else None
        model=read(sf/'active_model.mat',['value','model']);np.testing.assert_array_equal(model['value'],s['value'])
        full=float(np.linalg.norm(s['rhs']+model['model'])/np.linalg.norm(rhs));np.testing.assert_allclose(full,step['completeSelectedModelRelativeDefect'],rtol=1e-10)
        checks.append(dict(step=step['step'],directions=directions,activeArithmetic=arithmetic,projectorCompleteSupport=span_checks,publicTemperatureRows=public,completeModelRelativeDefect=full))
        del s,context,A,plan,data;gc.collect()
    return dict(steps=checks,allNativeBidiagonalizationColumnsIndependentlyVerified=False)


def main(mode,folder):
    assert mode in ('unit','real','inflight');folder=Path(folder).resolve();bindings={}
    def read(path,variables=None):
        path=Path(path);bindings[str(path.relative_to(ROOT))]=sha(path);return loadmat(path,simplify_cells=True,variable_names=variables)
    snapshot=None
    if mode=='inflight':
        progress_file=folder/'progress.json';digest=sha(progress_file)
        steps=json.loads(progress_file.read_text('utf8'));steps=[steps] if isinstance(steps,dict) else steps
        assert steps and sha(progress_file)==digest
        snapshot=folder/f'active_continuation_snapshot_step{len(steps):02d}.json';assert not snapshot.exists()
        snapshot.write_text(json.dumps(dict(status='frozen-inflight-completed-steps-only',steps=steps,
            originalProgressSHA256=digest,wholeRunTerminalClaimed=False),indent=2)+'\n',encoding='utf8')
        bindings[str(snapshot.relative_to(ROOT))]=sha(snapshot)
    checks=unit(folder,read) if mode=='unit' else real(folder,read,snapshot)
    assert all(sha(ROOT/path)==h for path,h in bindings.items())
    out=dict(status='passed-independent-temperature-band-active-cycles-and-complete-objective',mode=mode,checks=checks,
        physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__),recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid())
    out['wholeRunTerminalClaimed']=mode!='inflight'
    target=folder/('independent_temperature_active.json' if mode!='inflight' else f'independent_temperature_active_step{len(steps):02d}.json')
    assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))


if __name__=='__main__':main(*sys.argv[1:3])
