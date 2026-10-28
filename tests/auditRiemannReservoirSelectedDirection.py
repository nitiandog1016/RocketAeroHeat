# -*- coding: utf-8 -*-
"""用公开压力曲线/Brent独立核对采样Jacobian及原NS方向证据。"""
import json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT,VAL,sha,exact,PUBLIC_BINDINGS

case=sys.argv[1];assert case.replace('_','').isalnum();target=VAL/(case+'_independent_reservoir_selected_direction.json');assert not target.exists()
file=VAL/case/'riemann_selected_direction_unit_evidence.mat';s=loadmat(file,simplify_cells=True);gas=s['gas'];g=gas['gamma'];R=gas['R'];cv=R/(g-1)
bindings={str(file.relative_to(ROOT)):sha(file),**PUBLIC_BINDINGS}
def primitive(q):
    rho=q[0];v=q[1:4]/rho;return np.r_[rho,v,(q[4]-.5*rho*(v@v))/(rho*cv)]

worst=0.;tested=0;excluded=0
for i,(ql,qr,n) in enumerate(zip(s['QL'],s['QR'],s['normal'])):
    pl=primitive(ql);pr=primitive(qr);base,_,region=exact(pl,pr,n,gas)
    assert region==s['sample']['sampleRegion'][i]
    for side in range(2):
        q=ql if side==0 else qr;p=pl if side==0 else pr
        speed=max(np.linalg.norm(p[1:4]),np.sqrt(g*R*p[4]));basis=np.r_[p[0],np.repeat(p[0]*speed,3),p[0]*speed**2]
        matrix=s['sample']['left' if side==0 else 'right'][:,:,i]
        for b in range(5):
            h=1e-6*basis[b];plus=q.copy();minus=q.copy();plus[b]+=h;minus[b]-=h
            vp,_,rp=exact(primitive(plus),pr,n,gas) if side==0 else exact(pl,primitive(plus),n,gas)
            vm,_,rm=exact(primitive(minus),pr,n,gas) if side==0 else exact(pl,primitive(minus),n,gas)
            if rp!=region or rm!=region:excluded+=1;continue
            difference=float(np.max(abs((vp-vm)/(2*h)-matrix[:,b]))/max(1,np.max(abs(matrix[:,b]))))
            worst=max(worst,difference);tested+=1
assert tested==150 and excluded==0 and worst<5e-7
old_file=VAL/'branchjac01/exact_branch_jacobian_evidence.mat';bindings[str(old_file.relative_to(ROOT))]=sha(old_file)
old=loadmat(old_file,simplify_cells=True,variable_names=['left','right'])
np.testing.assert_array_equal(s['left'],old['left']);np.testing.assert_array_equal(s['right'],old['right'])
cache=s['cache'];ids=np.asarray(cache['outsideFaceIDs'],int)-1;owner=np.asarray(cache['outsideOwner'],int)-1;boundary_worst=0.;boundary_columns=0
for j,(face_id,cell_id) in enumerate(zip(ids,owner)):
    p=primitive(s['Q'][cell_id]);n=s['cut']['openFaces']['normal'][face_id]
    axis=int(s['cut']['openFaces']['axis'][face_id])-1;names=['xmin','xmax','ymin','ymax','zmin','zmax'];b=s['boundary'][names[2*axis+1]]
    r=b['primitive'];_,_,region=exact(p,r,n,gas);basis=np.r_[p[0],np.repeat(np.sqrt(g*R*p[4]),3),p[4]]
    for column in range(5):
        h=1e-6*basis[column];plus=p.copy();minus=p.copy();plus[column]+=h;minus[column]-=h
        vp,_,rp=exact(plus,r,n,gas);vm,_,rm=exact(minus,r,n,gas);assert rp==rm==region
        expected=cache['boundaryJacobian'][:,column,j]
        error=float(np.max(abs((vp-vm)/(2*h)-expected))/max(1,np.max(abs(expected))));boundary_worst=max(boundary_worst,error);boundary_columns+=1
assert boundary_columns==40 and boundary_worst<5e-7
differences=[]
for value_key,evidence_key in [('fd','fdEvidence'),('half','halfEvidence')]:
    e=s[evidence_key];rebuilt=(e['ratePlus']-e['rateMinusOrBase'])/s['scale']/e['actualScalarStep']
    np.testing.assert_array_equal(rebuilt,s[value_key]);assert e['actualOriginalResidualEvaluations']==1
    differences.append(float(np.linalg.norm(s['value']-rebuilt)/max(1,np.linalg.norm(s['value']))))
assert max(differences)<2e-6 and np.linalg.norm(s['value']-s['uncached'])/max(1,np.linalg.norm(s['value']))<1e-12
assert all(sha(ROOT/p)==digest for p,digest in bindings.items())
out=dict(status='passed-independent-150-public-sample-Jacobian-columns-40-boundary-columns-and-original-NS-arithmetic',
    processID=os.getpid(),sampleColumns=tested,crossBranchColumnsExcluded=excluded,maximumPublicSampleFDDifference=worst,
    boundaryPrimitiveColumns=boundary_columns,maximumPublicBoundaryPrimitiveFDDifference=boundary_worst,
    originalNSDirectionalRelativeDifferences=differences,oldFluxJacobianArraysExactlyPreserved=True,
    fullIndependentNSJacobianClaimed=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
