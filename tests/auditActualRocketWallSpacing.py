"""辨别层流原生Y+占位与基于实际首层距离/壁剪的网格准备估计。"""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v119_validation'
sys.path.insert(0,str(ROOT/'tests'));from auditSU2V110 import vtk_tokens

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'wall_spacing.json';assert not target.exists()
    mesh=ROOT/'outputs/v118_validation/actual01/mesh/hybrid_volume.mat'
    d=loadmat(str(mesh),squeeze_me=True,struct_as_record=False);V=d['sourceVertices'];F=d['sourceFaces'].astype(int)-1;n=len(V)
    nodes=d['vertices'];assert np.array_equal(nodes[:n],V)
    distance=np.linalg.norm(nodes[n:2*n]-V,axis=1)
    path=OUT/'continuation01/surface_flow.vtk';x,f=vtk_tokens(path);err,ids=cKDTree(V).query(x)
    assert max(err)<1e-6 and len(set(ids))==n
    cfg=json.loads((OUT/'continuation01/case_manifest.json').read_text())['flight']
    rhoinf=cfg['pressure']/(287.05287*cfg['temperature']);speed=cfg['Mach']*np.sqrt(1.4*287.05287*cfg['temperature'])
    qinf=.5*rhoinf*speed**2
    shear=np.linalg.norm(f['Skin_Friction_Coefficient'],axis=1)*qinf
    rho=f['Density'][:,0];mu=f['Laminar_Viscosity'][:,0]
    assert min(rho)>0 and min(mu)>0
    estimate=distance[ids]*np.sqrt(rho*shear)/mu
    node=np.zeros(n);node[ids]=estimate;face=node[F].mean(axis=1)
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    geometry=ROOT/'tools/runtime/su2-thermal-build/source/Common/src/geometry/CGeometry.cpp'
    forces=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl'
    assert np.all(f['Y_Plus']==0)
    result=dict(physicsQualified=False,newQualifiedArea=0,reportedYPlusAllZero=True,
                reason='fixed source computes nearest-neighbor distance for RANS-family solvers; laminar NS skips it, so native zero is not a mesh certificate',
                formula='h_first * sqrt(rho_wall * norm(Cf)*q_inf) / mu_wall; laminar saved shear only',
                firstLayerDistanceRangeM=[float(min(distance)),float(max(distance))],dynamicPressurePa=float(qinf),
                estimatedWallUnitsPercentiles=dict(zip(['min','p50','p95','p99','max'],np.percentile(estimate,[0,50,95,99,100]).tolist())),
                faceMeanEstimateAbove1AreaFraction=float(sum(area[face>1])/sum(area)),
                limitation='explicit first-normal-point estimate, not native RANS y+, not converged turbulent shear or heat-gradient validation',
                sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'tests/auditSU2V110.py',mesh,path,geometry,forces]})
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
