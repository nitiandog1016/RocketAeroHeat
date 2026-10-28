# -*- coding: utf-8 -*-
"""分解原NS实际失败方向的剩余量，并测量既有粗像还能去除多少，不新增NS调用。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/case/'coarse_direction_diagnostic.json';assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'wall_velocity_krylov_evidence.mat';bindings={str(file.relative_to(ROOT)):sha(file)}
s=loadmat(file,simplify_cells=True,variable_names=['result']);r=s['result'];logs=r['linearLog'];a=logs if isinstance(logs,dict) else list(logs)[-1]
e=a['steadyDirectionEvidence'];rhs=r['initialRate']/r['residualScale'];remaining=rhs+e['value'];reference=float(np.linalg.norm(rhs));norm=float(np.linalg.norm(remaining))
assert abs(norm/reference-a['fullSystemRelativeDefect'])<3e-12
source=Path(r['options']['coarseOperatorFile']);source=source if source.is_absolute() else ROOT/source;bindings[str(source.relative_to(ROOT))]=sha(source)
coarse=loadmat(source,simplify_cells=True,variable_names=['coarse'])['coarse'];T=coarse['image'];normal_rhs=T.T@remaining.ravel(order='F')
c=np.linalg.solve(coarse['normalMatrix'],normal_rhs);projected=remaining.ravel(order='F')-T@c
before2=norm**2;after2=float(np.dot(projected,projected));fraction=max(0,1-after2/before2)
nf=len(coarse['space']['geometry']['inputGeometry']['F']);nl=len(rhs)//nf;layer=remaining.reshape((nl,nf,5))
layer_share=np.sum(layer**2,axis=(1,2))/before2;equation_share=np.sum(remaining**2,axis=0)/before2;face_energy=np.sum(layer**2,axis=(0,2));order=np.argsort(face_energy)[-12:][::-1]
geometry=coarse['space']['geometry']['inputGeometry'];vertices=geometry['V'];faces=geometry['F'].astype(int)-1;centroids=vertices[faces].mean(axis=1)
normals=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]]);normals/=np.linalg.norm(normals,axis=1)[:,None]
axial_fraction=(centroids[:,0]-vertices[:,0].min())/(vertices[:,0].max()-vertices[:,0].min());axial_bin=np.minimum(9,np.floor(10*axial_fraction).astype(int))
axial_share=np.bincount(axial_bin,weights=face_energy,minlength=10)/before2
coarse_overlap=float(np.linalg.norm(normal_rhs)/max(1,norm));history=np.atleast_1d(a['gmresResidualHistory'])/reference
out=dict(status='source-bound-original-NS-failed-direction-residual-and-unused-coarse-projection-diagnostic',processID=os.getpid(),
    actualNSRelativeDefect=norm/reference,modelFinalRelativeDefect=float(history[-1]),modelActualRelativeDefectDifference=abs(float(history[-1])-norm/reference),
    coarseProjectionRemainingRelativeDefect=float(np.linalg.norm(projected)/reference),remainingEnergyRemovableByExistingCoarseFraction=fraction,
    normalizedCoarseImageInnerProducts=coarse_overlap,layerRemainingSquaredNormFractions=layer_share.tolist(),
    equationNames=['density','xMomentum','yMomentum','zMomentum','totalEnergy'],equationRemainingSquaredNormFractions=equation_share.tolist(),
    axialBinEdgesFromAftToNose=np.linspace(0,1,11).tolist(),axialRemainingSquaredNormFractions=axial_share.tolist(),
    dominantOriginalFaces=[dict(originalFaceID=int(i+1),centroidM=centroids[i].tolist(),bodyOutwardNormal=normals[i].tolist(),remainingSquaredNormFraction=float(face_energy[i]/before2)) for i in order],
    newActualNSEvaluations=0,newModelApplications=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
assert all(sha(ROOT/p)==h for p,h in bindings.items());target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
