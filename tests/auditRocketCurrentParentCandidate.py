# -*- coding: utf-8 -*-
"""独立核对当前C的LSQR坐标提升、完整像、原NS新方向和资格分离。"""
from pathlib import Path
import json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketParentImageScale import setup
from screenRocketFullFLeastSquares import ROOT,VAL,sha,local

case=sys.argv[1];assert case.replace('_','').isalnum();folder=VAL/case;target=VAL/(case+'_independent_current_candidate.json');assert not target.exists();bindings={}
def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)
data=read(folder/'current_parent_candidate_evidence.mat');packed=read(data['packedFile'],['packed'])['packed'];factor=read(data['factorFile'],['factor'])['factor']
A,plan,order=setup(packed,factor);del packed,factor
scale=read(data['scaleFile'],['columnScale'])['columnScale'];n=len(order)
correction=local(plan,(data['coordinates'][order]*scale[order]).ravel()).reshape(n,5);lift=data['initial'].copy();lift[order]+=correction
lift_error=float(np.linalg.norm(lift-data['value'])/max(1,np.linalg.norm(data['value'])));assert lift_error<1e-8
image=np.empty_like(lift);image[order]=(A@data['value'][order].ravel()).reshape(n,5)
image_error=float(np.linalg.norm(image-data['image'])/max(1,np.linalg.norm(data['image'])));assert image_error<1e-12
residual=float(np.linalg.norm(data['rhs']-image)/np.linalg.norm(data['rhs']));assert abs(residual-data['detail']['relativeCResidual'])<1e-12
assert data['detail']['internalCTargetMet']==(residual<=1e-3) and data['detail']['options']['relativeTolerance']==1e-3
assert data['detail']['nativeIterations']==256
model=read(folder/'current_parent_model_probe.mat');np.testing.assert_array_equal(model['value'],data['value'])
full=float(np.linalg.norm(data['rhs']+model['fullModel'])/np.linalg.norm(data['rhs']));assert abs(full-model['report']['fullSelectedModelOriginalRHSRelativeResidual'])<1e-12
probe=read(VAL/'directnsprobereal03/direct_ns_direction_probe.mat');np.testing.assert_array_equal(probe['directionU'],data['value'])
rhs=probe['baseRate']/probe['scale'];np.testing.assert_array_equal(rhs,data['rhs']);e=probe['nsEvidence']
jv=(e['ratePlus']-e['rateMinusOrBase'])/probe['scale']/e['actualScalarStep'];np.testing.assert_array_equal(jv,e['value'])
actual=float(np.linalg.norm(rhs+jv)/np.linalg.norm(rhs));assert actual<=.1 and abs(actual-probe['report']['originalNSDirectionDefect'])<1e-12
source=read(data['sourceBinding']['sourceFile'],['result'])['result'];np.testing.assert_array_equal(probe['Q'],source['Q'])
assert not model['report']['originalNSDirectionCheckedHere'] and not model['report']['physicalStateChanged']
assert probe['report']['newFlowAcceptances']==0 and model['report']['newProductionSolvedAreaM2']==0
assert not data['detail']['internalCTargetMet']
assert all(sha(ROOT/p)==digest for p,digest in bindings.items())
out=dict(status='passed-independent-current-parent-LSQR-coordinate-lift-full-image-and-separate-original-NS-direction',
    processID=os.getpid(),originalCVs=n,nativeIterations=256,coordinateLiftRelativeDifference=lift_error,fullCImageRelativeDifference=image_error,
    currentCRelativeResidual=residual,currentCTargetStillNotMet=True,fullSelectedModelRelativeResidual=full,
    separateOriginalNSDirectionDefect=actual,separateOriginalNSPointOneMet=True,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
