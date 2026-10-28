# -*- coding: utf-8 -*-
"""核对真实Q20全部储库采样Jacobian及完整NS方向差分，不声称独立NS求解。"""
import json,os,sys,math
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT,VAL,sha,exact,PUBLIC_BINDINGS

case=sys.argv[1];assert case.replace('_','').isalnum();target=VAL/(case+'_independent_real_reservoir_direction.json');assert not target.exists()
file=VAL/case/'riemann_selected_direction_real_evidence.mat';s=loadmat(file,simplify_cells=True);gas=s['gas'];scale=s['scale'];rhs=s['rate']/scale
bindings={str(file.relative_to(ROOT)):sha(file),**PUBLIC_BINDINGS};source=s['report']['sourceFile'];assert sha(source)==s['report']['sourceSHA256']
bindings[str(__import__('pathlib').Path(source).relative_to(ROOT))]=sha(source)
differences=[];defects=[]
values=s['fdValues'];evidence=s['fdEvidence']
for value,e in zip(values,evidence):
    np.testing.assert_array_equal(e['directionU'],s['direction']);np.testing.assert_array_equal(e['rateMinusOrBase'],s['rate'])
    rebuilt=(e['ratePlus']-e['rateMinusOrBase'])/scale/e['actualScalarStep'];np.testing.assert_array_equal(rebuilt,value)
    differences.append(float(np.linalg.norm(s['value']-rebuilt)/max(1,np.linalg.norm(s['value']))))
    defects.append(float(np.linalg.norm(rhs+rebuilt)/np.linalg.norm(rhs)))
np.testing.assert_allclose(differences,s['report']['analyticToOriginalNSRelativeDifferences'],rtol=3e-13,atol=1e-15)
accurate_norm=lambda v:math.sqrt(math.fsum(np.asarray(v).ravel()**2))
accurate_defects=[accurate_norm(rhs+v)/accurate_norm(rhs) for v in values]
# 百万项正平方和的范数实现舍入界；独立参考采用fsum，原NS0.1门槛不变。
norm_roundoff_bound=4*rhs.size*np.finfo(float).eps/(1-rhs.size*np.finfo(float).eps)
norm_reference_difference=float(np.max(abs(np.asarray(accurate_defects)-s['report']['originalNSRelativeDefects'])/np.maximum(1,abs(s['report']['originalNSRelativeDefects']))))
assert norm_reference_difference<=norm_roundoff_bound
assert max(differences)<2e-5 and min(defects)>.1 and s['report']['newFlowAcceptances']==0
cache=s['cache'];geometry=cache['geometry'];ids=np.asarray(cache['outsideFaceIDs'],int)-1;owners=np.asarray(cache['outsideOwner'],int)-1
P=cache['state']['interiorPrimitive'];matrices=cache['boundaryJacobian'];saved=cache['state']['primitive'][ids];incoming=cache['state']['incoming'][ids]
del cache;s.pop('cache')
counts=0;excluded=0;worst=0.;base_error=0.;worst_location=None
for j,(face_id,owner) in enumerate(zip(ids,owners)):
    p=P[owner];normal=geometry['openFaces']['normal'][face_id];axis=int(geometry['openFaces']['axis'][face_id])-1
    name=['xmax','ymax','zmax'][axis];reservoir=s['boundary'][name]['primitive'];base,_,region=exact(p,reservoir,normal,gas)
    base_error=max(base_error,float(np.max(abs(base-saved[j]))/max(1,np.max(abs(base)))))
    sound=np.sqrt(gas['gamma']*gas['R']*p[4]);out_sound=np.sqrt(gas['gamma']*gas['R']*base[4])
    input_scale=np.r_[p[0],np.repeat(sound,3),p[4]];output_scale=np.r_[base[0],np.repeat(out_sound,3),base[4]]
    for column in range(5):
        h=1e-6*input_scale[column];plus=p.copy();minus=p.copy();plus[column]+=h;minus[column]-=h
        vp,_,rp=exact(plus,reservoir,normal,gas);vm,_,rm=exact(minus,reservoir,normal,gas)
        if rp!=region or rm!=region:
            excluded+=1;continue
        actual=(vp-vm)/(2e-6*output_scale);expected=matrices[:,column,j]*input_scale[column]/output_scale
        error=float(np.max(abs(actual-expected))/max(1,np.max(abs(expected))))
        if error>worst:worst=error;worst_location=dict(openFaceID=int(face_id+1),ownerCellID=int(owner+1),column=column+1)
        counts+=1
assert counts+excluded==len(ids)*5 and counts>98000 and worst<2e-6 and base_error<3e-12
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-real-Q20-full-NS-step-arithmetic-and-all-reservoir-boundary-columns',processID=os.getpid(),
    boundaryFaces=len(ids),dimensionlessStableBranchColumns=counts,crossBranchColumnsExcluded=excluded,
    maximumScaledPublicBoundaryFDDifference=worst,worstColumnLocation=worst_location,maximumPublicBaseTraceDifference=base_error,
    derivativeScaleDefinition='input/output density, sound-speed for three velocities, and temperature; epsilon=1e-6',
    accurateNSDefects=accurate_defects,normSummationReferenceRelativeDifference=norm_reference_difference,normSummationForwardBound=norm_roundoff_bound,
    originalNSDefects=defects,analyticToOriginalNSDifferences=differences,directionPointOneTargetStillNotMet=True,
    fullIndependentNSJacobianClaimed=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
