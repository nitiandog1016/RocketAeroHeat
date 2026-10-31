"""按SU2中位对偶定义重构节点体积，分开比较积分残差及局部变化率。"""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from diagnoseSU2Checkpoint import read_volume

ROOT=Path(__file__).resolve().parents[1]
FACES={12:[[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]],
       13:[[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]]}


def dual_volumes(v):
    volume=np.zeros(len(v.vertices));cell_sum=np.zeros(len(v.cells));minimum=np.inf
    for kind,faces in FACES.items():
        selected=np.flatnonzero(v.cellTypes==kind)
        if not len(selected):continue
        ids=np.stack(v.cells[selected]).astype(int)-1
        q=v.vertices[ids];center=q.mean(axis=1)
        for face in faces:
            fc=q[:,face].mean(axis=1)
            for a,b in zip(face,face[1:]+face[:1]):
                midpoint=(q[:,a]+q[:,b])/2
                for k in (a,b):
                    p=q[:,k];tet=np.abs(np.einsum('ij,ij->i',np.cross(midpoint-p,fc-p),center-p))/6
                    assert np.all(tet>0);minimum=min(minimum,float(tet.min()))
                    np.add.at(volume,ids[:,k],tet);cell_sum[selected]+=tet
    assert np.all(volume>0)
    error=float(np.max(abs(cell_sum/v.cellVolume-1)))
    assert error<2e-10 and abs(volume.sum()/v.fluidVolume-1)<2e-12
    return volume,dict(maxCellVolumeRelativeDifference=error,domainVolume=float(volume.sum()),minimumSubVolume=minimum)


def metrics(xyz,fields,volume,wall_count):
    r=np.column_stack((fields['Residual_Density'],fields['Residual_Momentum'],fields['Residual_Energy']))
    masks={'wall':np.arange(len(xyz))<wall_count,'normalShell':(np.arange(len(xyz))//wall_count>0)&(np.arange(len(xyz))//wall_count<=24),
           'outer':np.arange(len(xyz))//wall_count>24,'downstream':xyz[:,0]<0}
    norms={'integratedSquare':r*r,'localRateSquare':(r/volume[:,None])**2,'volumeWeightedRateSquare':r*r/volume[:,None]}
    result={}
    for name,square in norms.items():
        total=square.sum(axis=0);result[name]={'fractions':{k:(square[m].sum(axis=0)/total).tolist() for k,m in masks.items()},
                                            'sum':total.tolist(),'maxEnergyNode':int(np.argmax(square[:,-1]))}
    return result


def main():
    path=ROOT/'outputs/v108_validation/normal01/normal_volume.mat';v=loadmat(str(path),squeeze_me=True,struct_as_record=False)['volume']
    volume,check=dual_volumes(v);source=ROOT/'outputs/v110_validation/probe02/volume_flow.vtk';xyz,f=read_volume(source)
    d,ids=cKDTree(v.vertices).query(xyz);assert max(d)<1e-6 and len(set(ids))==len(ids)
    order=np.argsort(ids);f={k:a[order] for k,a in f.items()};nw=len(np.unique(np.concatenate(v.wall)))
    report=dict(geometry=check,normalization=metrics(v.vertices,f,volume,nw),
                semantics='R integrated flux balance; R/V local rate; sum(R^2/V) volume-weighted squared rate; none is physical error',
                inputSHA256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in (path,source,Path(__file__).resolve())})
    out=ROOT/'outputs/v111_validation';out.mkdir(parents=True,exist_ok=True)
    assert not (out/'dual_volume.npy').exists();np.save(str(out/'dual_volume.npy'),volume)
    (out/'normalization.json').write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report,indent=2))


if __name__=='__main__':main()
