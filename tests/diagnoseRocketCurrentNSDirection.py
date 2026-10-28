# -*- coding: utf-8 -*-
"""比较两个真实基态对同一方向的原NS作用；仅定位变化，不提出接受或面积资格。"""
from pathlib import Path
import hashlib,json,os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'current_Q20_NS_direction_diagnostic.json';assert not target.exists()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
bindings={}
def read(path,variables):
    bindings[path.relative_to(ROOT).as_posix()]=sha(path);return loadmat(path,simplify_cells=True,variable_names=variables)
a=read(VAL/'directnsprobereal01/direct_ns_direction_probe.mat',['Q','baseRate','scale','directionU','nsEvidence'])
b=read(VAL/'directnsprobereal02/direct_ns_direction_probe.mat',['Q','baseRate','scale','directionU','nsEvidence'])
np.testing.assert_array_equal(a['directionU'],b['directionU']);np.testing.assert_array_equal(a['scale'],b['scale'])
np.testing.assert_array_equal(b['Q'],a['Q']+(1/256)*a['directionU']*a['scale'])
domain=read(ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat',['domain'])['domain'];nf=len(domain['inputGeometry']['F'])
delta=b['nsEvidence']['value']-a['nsEvidence']['value'];rhs=b['baseRate']/b['scale'];norm=float(np.linalg.norm(rhs))
energy=np.sum(delta**2,axis=1);order=np.argsort(energy)[::-1];total=float(energy.sum());nl=len(energy)//nf
details=[dict(cellID=int(i+1),layer=int(i//nf+1),sourceFacetID=int(i%nf+1),centroidM=domain['cells']['centroid'][i].tolist(),
              squaredChangeFraction=float(energy[i]/total)) for i in order[:20]]
old=read(VAL/'riemannreservoirreal01/riemann_reservoir_real_evidence.mat',['info'])['info']['riemannReservoirBoundaryAudit']
new=read(VAL/'krylov20/wall_velocity_krylov_evidence.mat',['result'])['result']['info']['riemannReservoirBoundaryAudit']
np.testing.assert_array_equal(old['faceID'],new['faceID'])
flow=np.flatnonzero(np.asarray(old['incoming'],bool)!=np.asarray(new['incoming'],bool));wave=np.flatnonzero(old['sampleRegion']!=new['sampleRegion'])
changed_faces=np.asarray(new['faceID'],int)[flow]-1
owners=np.asarray(domain['openFaces']['cellMinus'],int)[changed_faces]-1
out=dict(status='diagnosed-same-direction-current-Q20-original-NS-variation',processID=os.getpid(),
    originalNSActionChangeRelativeToCurrentRHS=float(np.linalg.norm(delta)/norm),
    layerSquaredChangeFractions=(energy.reshape(nl,nf).sum(axis=1)/total).tolist(),
    equationSquaredChangeFractions=(np.sum(delta**2,axis=0)/total).tolist(),topCells=details,
    changedBoundaryMassDirections=len(flow),changedBoundaryWaveRegions=len(wave),
    changedMassDirectionFaces=[dict(openFaceID=int(face+1),ownerCellID=int(owner+1),centroidM=domain['openFaces']['centroid'][face].tolist()) for face,owner in zip(changed_faces,owners)],
    first20CellsSquaredChangeFraction=float(energy[order[:20]].sum()/total),
    sourceDefectBefore=float(np.linalg.norm(a['baseRate']/a['scale']+a['nsEvidence']['value'])/np.linalg.norm(a['baseRate']/a['scale'])),
    sourceDefectCurrent=float(np.linalg.norm(rhs+b['nsEvidence']['value'])/norm),
    physicalStateAcceptedHere=False,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','topCells','changedMassDirectionFaces')}))
