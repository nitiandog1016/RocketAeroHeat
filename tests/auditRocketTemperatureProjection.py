"""独立核对公开Riemann温度导数、稠密受限解、完整坐标提升和原目标。"""
from pathlib import Path
from datetime import datetime,timezone
import os,sys,json,gc
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import numpy as np
from scipy.io import loadmat
from scipy.linalg import null_space
from screenRocketFullFLeastSquares import ROOT,sha,local
from screenRocketParentImageScale import setup
from screenRocketAugmentedPatchQR import metric_blocks,block_weight
from riemannReservoirReference import exact,PUBLIC_BINDINGS


def row_reference(rows,step_factor=1e-6):
    # 可显式细化差分步长以核对高动能/低内能制造场；缺省及1e-7核对门槛保持。
    assert np.isfinite(step_factor) and step_factor>0
    gas=rows['gas'];ql=np.atleast_2d(rows['QL']);qr=np.atleast_2d(rows['QR']);normals=np.atleast_2d(rows['normal'])
    scale=rows['scale'];analytic=np.column_stack([np.atleast_2d(rows['leftRowsK']),np.atleast_2d(rows['rightRowsK'])])
    fd=np.zeros_like(analytic);base=[]
    def primitive(q):
        rho=q[0];u=q[1:4]/rho;T=(q[4]-.5*rho*(u@u))/(rho*gas['R']/(gas['gamma']-1))
        return np.r_[rho,u,T]
    for i,(l,r,n) in enumerate(zip(ql,qr,normals)):
        s,_,region=exact(primitive(l),primitive(r),n,gas);base.append(s[4])
        assert region==np.atleast_1d(rows['sampleRegion'])[i]
        for side in range(2):
            for b in range(5):
                pair=np.row_stack([l,r]);h=step_factor*max(1.,abs(pair[side,b]/scale[b]))
                plus=pair.copy();minus=pair.copy();plus[side,b]+=h*scale[b];minus[side,b]-=h*scale[b]
                sp,_,rp=exact(primitive(plus[0]),primitive(plus[1]),n,gas)
                sm,_,rm=exact(primitive(minus[0]),primitive(minus[1]),n,gas)
                assert rp==rm==region
                fd[i,side*5+b]=(sp[4]-sm[4])/(2*h)
    difference=float(np.linalg.norm(fd-analytic)/max(1.,np.linalg.norm(analytic)))
    assert difference<1e-7
    np.testing.assert_allclose(base,np.atleast_1d(rows['baseTemperatureK']),rtol=1e-12)
    return dict(allTemperatureDerivativeColumnCount=fd.size,publicTemperatureFDRelativeDifference=difference,
                baseTemperatureK=base,allPerturbationsStayedInSameSelectedRegion=True)


def unit(folder,read):
    s=read(folder/'temperature_projection_unit.mat',None);reference=row_reference(s['rows']);checks=[]
    for c in s['cases']:
        H=c['rows']['H'].toarray();b=c['source'];rhs=b['base']['rhs'].ravel(order='F');N=len(rhs)
        penalty=c['factor']['penaltyInOriginalResidualUnits'];A=np.row_stack([b['base']['C'],penalty*b['W']]);Z=null_space(H)
        dense=Z@np.linalg.lstsq(A@Z,np.r_[rhs,np.zeros(N)],rcond=None)[0]
        value=c['value'].ravel(order='F');error=float(np.linalg.norm(value-dense)/max(1.,np.linalg.norm(dense)));assert error<1e-8
        Pi=c['Pi'];P=c['P'];lift=P@Pi@c['coordinates'].ravel(order='F')
        assert np.linalg.norm(lift-value)/max(1.,np.linalg.norm(value))<1e-12
        assert np.linalg.norm(Pi-Pi.T)<1e-12 and np.linalg.norm(Pi@Pi-Pi)<1e-12
        residual=float(np.linalg.norm(H@value)/max(1.,np.linalg.norm(H)*np.linalg.norm(value)));assert residual<=1e-10
        checks.append(dict(denseConstrainedDifference=error,constraintNormalizedError=residual))
    return dict(publicDerivative=reference,cases=checks)


def real(folder,read):
    s=read(folder/'temperature_projection_directions.mat',None)
    f=read(folder/'temperature_projection_system.mat',['factor','projector'])
    pfile=Path(s['preparationFile']);assert sha(pfile)==s['preparationSHA256'];p=json.loads(pfile.read_text(encoding='utf-8-sig'))
    assert sha(p['sourceBinding']['sourceFile'])==p['sourceBinding']['sourceFileSHA256']
    for a in p['artifacts'].values():assert sha(a['file'])==a['sha256']
    packed=read(p['artifacts']['packed']['file'],['packed'])['packed'];factor=f['factor']
    buckets=[factor['buckets']] if isinstance(factor['buckets'],dict) else list(np.atleast_1d(factor['buckets']))
    for b in buckets:b['inverse']=b['map']
    factor['buckets']=buckets;A,plan,order=setup(packed,factor);del packed,factor,buckets
    gc.collect();n=len(order);H=s['rows']['H'];basis=np.asarray(f['projector']['basis']).reshape(5*n,-1,order='F')
    basis_algebraic=np.column_stack([basis[:,j].reshape(n,5,order='F')[order].ravel() for j in range(basis.shape[1])])
    np.testing.assert_allclose(basis.T@basis,np.eye(basis.shape[1]),rtol=0,atol=1e-12)
    expected_basis=[]
    for row in H.toarray():
        raw=row.reshape(n,5,order='F')[order].ravel();expected_basis.append(local(plan,raw,True))
    original=np.column_stack(expected_basis)
    relation=np.linalg.norm(original-basis_algebraic@(basis_algebraic.T@original))/np.linalg.norm(original);assert relation<1e-12
    W=metric_blocks(s['metric'],order);rhs=s['rhs'][order].ravel();penalty=f['factor']['penaltyInOriginalResidualUnits']
    results=[]
    for label in ('free','projected'):
        value=s[label][order].ravel();coordinates=s[label+'Coordinates'][order].ravel()
        if label=='projected':coordinates=coordinates-basis_algebraic@(basis_algebraic.T@coordinates)
        lifted=local(plan,coordinates);error=np.linalg.norm(lifted-value)/max(1.,np.linalg.norm(value));assert error<1e-12
        image=A@value;weighted=block_weight(W,value);relative=float(np.linalg.norm(rhs-image)/np.linalg.norm(rhs))
        objective=relative**2+(penalty*np.linalg.norm(weighted)/np.linalg.norm(rhs))**2
        detail=s[label+'Detail'];np.testing.assert_allclose(relative,detail['relativeCResidual'],rtol=1e-10)
        np.testing.assert_allclose(objective,detail['normalizedAugmentedObjective'],rtol=1e-10)
        derivative=np.atleast_1d(H@s[label].ravel(order='F'))
        if label=='projected':
            bound=np.linalg.norm(derivative)/max(1.,np.linalg.norm(H.data)*np.linalg.norm(value));assert bound<=1e-10
            np.testing.assert_allclose(derivative,np.atleast_1d(detail['constraintDerivativeK']),atol=1e-8,rtol=1e-8)
        branch=read(folder/label/'direction.mat',['value','model']);np.testing.assert_array_equal(branch['value'],s[label])
        model=float(np.linalg.norm(s['rhs']+branch['model'])/np.linalg.norm(rhs))
        results.append(dict(method=label,coordinateLiftRelativeDifference=float(error),originalCRelativeResidual=relative,
                            normalizedAugmentedObjective=float(objective),constraintDerivativeK=derivative.tolist(),
                            completeSelectedModelRelativeDefect=model))
    return dict(publicDerivative=row_reference(s['rows']),constraintRightCoordinateSpanRelativeDifference=float(relation),
                cases=results,allNativeBidiagonalizationColumnsIndependentlyVerified=False)


def main(mode,folder):
    assert mode in ('unit','real');folder=Path(folder).resolve();bindings={**PUBLIC_BINDINGS}
    def read(path,variables):
        path=Path(path);bindings[str(path.relative_to(ROOT))]=sha(path)
        return loadmat(path,simplify_cells=True,variable_names=variables)
    checks=unit(folder,read) if mode=='unit' else real(folder,read)
    assert all(sha(ROOT/path)==h for path,h in bindings.items())
    out=dict(status='passed-independent-public-temperature-derivative-projection-and-complete-objective',mode=mode,checks=checks,
             physicalSolution=False,newProductionSolvedAreaM2=0,scriptSHA256=sha(__file__),inputSHA256=bindings,
             recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid())
    target=folder/'independent_temperature_projection.json';assert not target.exists()
    target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))


if __name__=='__main__':main(*sys.argv[1:3])
