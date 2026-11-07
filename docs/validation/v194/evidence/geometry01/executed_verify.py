# -*- coding: utf-8 -*-
"""独立面积、体积、矩阵范数质量与局部约束复核；失败网格不授予运行资格。"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import time
import numpy as np
from scipy.io import loadmat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from readEilmerNative import sha
from tetgenFixedInterface import sorted_cells
from auditSU2HybridTopology import paired
from cgalTetrahedralContract import faces_and_neighbors
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parents[1]


def verify(case,out):
    """模块一：真值来自原壁和代数独立表达式，不用生成器质量函数作真值。"""
    start=time.monotonic();case,out=map(Path,(case,out));out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_verify.py')
    report=json.loads((case/'report.json').read_text())
    for n,d in report['filesSHA256'].items():assert sha(case/n)==d
    m=np.load(case/'mesh.npz');P=m['vertices'];T=m['tetrahedra'];F=m['wall'];G=m['farfield']
    original=loadmat(ROOT/'outputs/v104_validation/baseline02/current_case.mat',simplify_cells=True)['current']
    assert np.array_equal(P[:9878],original['V']) and np.array_equal(F+1,original['F'])
    assert np.array_equal(m['originalFaceIDs'],np.arange(19752)+1)
    q=P[F];lengths=np.linalg.norm(q-np.roll(q,1,axis=1),axis=2);s=lengths.sum(1)/2
    area=np.sqrt(s*np.prod(s[:,None]-lengths,axis=1))
    assert abs(area.sum()-1.9646292565715551)<2e-13
    # 用行列式计算体积，用边界散度积分核对外壳减内壳。
    tet=P[T];A=tet[:,1:]-tet[:,:1];volume=np.linalg.det(A)/6
    assert np.all(volume>0)
    boundary=np.r_[F[:,[0,2,1]],G]
    poly=float(np.linalg.det(P[boundary]).sum()/6)
    assert abs(volume.sum()/poly-1)<2e-10
    actual,graph=faces_and_neighbors(T)
    assert np.array_equal(sorted_cells(actual),sorted_cells(boundary))
    pairing=paired(T[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3),boundary)
    assert connected_components(graph,directed=False,return_labels=False)==1
    # q=3|det J|/(||J||F ||cofactor J||F)，独立于生成器SVD表达式。
    equilateral=np.array([[1,.5,.5],[0,np.sqrt(3)/2,np.sqrt(3)/6],[0,0,np.sqrt(2/3)]])
    J=A.transpose(0,2,1)@np.linalg.inv(equilateral)
    cof=np.stack([np.cross(J[:,:,1],J[:,:,2]),np.cross(J[:,:,2],J[:,:,0]),np.cross(J[:,:,0],J[:,:,1])],axis=2)
    quality=3*np.linalg.det(J)/np.sqrt(np.sum(J*J,axis=(1,2))*np.sum(cof*cof,axis=(1,2)))
    assert abs(quality.min()-report['minimumSICN'])<1e-12
    bad=np.flatnonzero(quality<.1);fixed=set(np.unique(boundary).tolist());wallkeys={tuple(sorted(f)):i for i,f in enumerate(F)}
    rows=[]
    for i in bad:
        # 每个固定壁面给出的可达上界只证明该单面约束，不证明全域可修复。
        constraints=[]
        for order in ([0,1,2],[0,1,3],[0,2,3],[1,2,3]):
            ids=T[i,order];key=tuple(sorted(ids))
            if key in wallkeys:
                tri=P[ids];a=tri[1]-tri[0];b=(tri[2]-tri[0]-a/2)*2/np.sqrt(3)
                gram=np.array([[a@a,a@b],[a@b,b@b]])
                lam=np.linalg.eigvalsh(gram);k=np.sqrt(lam[1]/lam[0])
                constraints.append(dict(originalFaceID=wallkeys[key]+1,isolatedFaceQualityUpperBound=float(3/(1+k+1/k))))
        rows.append(dict(cellOneBased=int(i+1),quality=float(quality[i]),centerM=tet[i].mean(0).tolist(),
                         vertexIDs=T[i].tolist(),fixedVertexCount=sum(int(n) in fixed for n in T[i]),wallConstraints=constraints))
    np.savetxt(out/'bad_cells.csv',np.array([[r['cellOneBased'],r['quality'],r['fixedVertexCount']]+r['centerM'] for r in rows]).reshape(-1,6),
               delimiter=',',comments='',header='cell_one_based,SICN,fixed_vertices,x_m,y_m,z_m',fmt='%.17g')
    result=dict(passed=True,scope='independent geometry verification, not candidate acceptance',originalWallExact=True,
        originalAreaM2=float(area.sum()),volumeM3=float(volume.sum()),boundaryVolumeM3=poly,
        minimumVolumeM3=float(volume.min()),minimumSICN=float(quality.min()),badCells=rows,
        geometryQualityAccepted=bool(len(bad)==0),pairing=pairing,physicalQualification=False,
        inputsSHA256={str(p):sha(p) for p in [case/'mesh.npz',case/'report.json',Path(__file__)]},
        elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['passed','geometryQualityAccepted','minimumSICN','badCells','elapsedSeconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();verify(a.case,a.out)
