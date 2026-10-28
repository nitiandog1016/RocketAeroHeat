# -*- coding: utf-8 -*-
"""完整父C同端点/32方向的固定像列归一GMRES对照，避免法方程条件数平方。"""
from pathlib import Path
import os,json
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat,savemat
from scipy.sparse.linalg import gmres
from screenRocketFullFLeastSquares import ROOT,VAL,sha,local
from screenRocketParentImageScale import setup,make_operator

folder=VAL/'parentCimagegmresscreen02';assert not folder.exists();folder.mkdir();bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


p=read(VAL/'coupledpatchreal01/coupled_patch_factor_evidence.mat',['factor','packedFile','sourceBinding'])
assert sha(p['packedFile'])=='beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'
packed=read(p['packedFile'],['packed'])['packed'];A,plan,order=setup(packed,p['factor'])
s=read(VAL/'parentCimagescalescreen02/parent_image_scale_evidence.mat',['columnScaleAlgebraic','originalOrder'])
np.testing.assert_array_equal(s['originalOrder'],order+1)
base=read(VAL/'spatialreductionreal01/spatial_reduction_real_evidence.mat',['rhs','value'])
b=base['rhs'][order].ravel();initial=base['value'][order].ravel();r0=b-A@initial
B,counts=make_operator(A,plan,s['columnScaleAlgebraic']);history=[]
coordinates,info=gmres(B,r0,tol=1e-3*np.linalg.norm(b)/np.linalg.norm(r0),atol=0.,restart=32,maxiter=1,callback=history.append,callback_type='pr_norm')
value=initial+local(plan,s['columnScaleAlgebraic'].ravel()*coordinates);image=A@value
residual=float(np.linalg.norm(b-image)/np.linalg.norm(b));physical=np.empty_like(base['rhs']);physical[order]=value.reshape(len(order),5)
physical_image=np.empty_like(physical);physical_image[order]=image.reshape(len(order),5)
reference_file=VAL/'parentCimagescalescreen02/parent_image_scale_screen_summary.json';reference=json.loads(reference_file.read_text())
bindings[reference_file.relative_to(ROOT).as_posix()]=sha(reference_file)
savemat(folder/'parent_image_gmres_screen_evidence.mat',dict(rhs=base['rhs'],value=physical,image=physical_image,sourceBinding=p['sourceBinding']),do_compression=True)
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='tested-same-complete-parent-C-endpoint-fixed-image-scaled-GMRES-32-comparison',processID=os.getpid(),actualArnoldiSteps=len(history),nativeInfo=int(info),
    initialOriginalCResidual=float(np.linalg.norm(r0)/np.linalg.norm(b)),finalOriginalCResidual=residual,originalCTargetMet=residual<=1e-3,
    sameEndpoint32StepScaledLSQRResidual=reference['comparisons'][1]['finalOriginalCResidual'],operatorApplications=counts,
    originalParentResidualNormUnchanged=True,oldSchurUsed=False,fixedCRepresentativesUsed=False,
    allOriginalCVs=len(order),allFiveConservativeFieldsRetained=True,physicalStateChanged=False,
    independentAllRealArnoldiColumnsClaimed=False,newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
(folder/'parent_image_gmres_screen_summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
