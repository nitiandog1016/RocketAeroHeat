"""把实际原区域棱柱连接转换为polyMesh；不移动节点、不增删近壁层。"""
import time,shutil,json
import numpy as np
from scipy.io import loadmat
from prepareOpenFoamControl import ROOT,sha,save,header

BASE=ROOT/'outputs/v208_validation'

def main():
    start=time.monotonic();out=BASE/'layer_mesh01';out.mkdir(exist_ok=False)
    source=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';selection=BASE/'layers02/selection.npz'
    assert sha(source)=='22b35b72396c1a7d14195af9bd24d3c6328138484d4009638c1ab48cddec1d6c'
    data=loadmat(str(source),variable_names=['vertices','wedges','sourceFaces']);s=np.load(selection)
    original_rows=s['originalWedgeRows'];W=data['wedges'][original_rows].astype(np.int64)-1
    old_points=np.unique(W);inverse=np.full(len(data['vertices']),-1,np.int64);inverse[old_points]=np.arange(len(old_points))
    P=data['vertices'][old_points];C=inverse[W];nc=len(C);nv=9878
    original_faces=s['originalWallFaceIDs'];wall=inverse[data['sourceFaces'][original_faces-1]-1]
    top=inverse[data['sourceFaces'][original_faces-1]-1+32*nv]
    # 模块一：按全部节点编号精确配对，保留四边形，不将薄棱柱拆为tet。
    order=[[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]]
    raw=np.full((5*nc,4),-1,np.int64);counts=np.repeat([3,3,4,4,4],nc)
    for i,idx in enumerate(order):raw[i*nc:(i+1)*nc,:len(idx)]=C[:,idx]
    keys=np.sort(raw,axis=1);unique,first,inv,multiplicity=np.unique(keys,axis=0,return_index=True,return_inverse=True,return_counts=True)
    assert multiplicity.max()==2
    raw_owner=np.tile(np.arange(nc),5);owner=np.full(len(unique),nc,np.int64);other=np.full(len(unique),-1,np.int64)
    np.minimum.at(owner,inv,raw_owner);np.maximum.at(other,inv,raw_owner)
    F=raw[first];sizes=counts[first];internal=multiplicity==2
    assert np.all(owner[internal]<other[internal])
    boundary=~internal
    # 连接键准确识别原壁与外层，其余边界只能是原切口的层列四边面。
    wall_keys={tuple(sorted([-1,*t])):i for i,t in enumerate(wall)}
    top_keys={tuple(sorted([-1,*t])) for t in top}
    wall_rows=[];wall_ids=[];top_rows=[];cut_rows=[]
    for i in np.flatnonzero(boundary):
        key=tuple(unique[i])
        if key in wall_keys:wall_rows.append(i);wall_ids.append(original_faces[wall_keys[key]])
        elif key in top_keys:top_rows.append(i)
        else:cut_rows.append(i)
    assert len(wall_rows)==len(top_rows)==7838 and len(cut_rows)==46*32
    assert np.all(sizes[cut_rows]==4)
    internal_rows=np.flatnonzero(internal)
    internal_rows=internal_rows[np.lexsort((other[internal_rows],owner[internal_rows]))]
    idx=np.r_[internal_rows,wall_rows,top_rows,cut_rows];ni=int(internal.sum())
    F=F[idx];sizes=sizes[idx];owner=owner[idx];neighbour=other[idx[:ni]]
    centers=P[C].mean(axis=1);sf=np.zeros((len(F),3));fc=np.zeros((len(F),3))
    for count in [3,4]:
        mask=sizes==count;q=P[F[mask,:count]]
        fc[mask]=q.mean(axis=1);sf[mask]=np.cross(q,np.roll(q,-1,axis=1)).sum(axis=1)/2
        flip=np.einsum('ij,ij->i',sf[mask],fc[mask]-centers[owner[mask]])<0
        rows=np.flatnonzero(mask)[flip];F[rows,:count]=F[rows,:count][:,::-1];sf[rows]*=-1
    assert np.all(np.einsum('ij,ij->i',sf,fc-centers[owner])>0)
    assert np.all(np.einsum('ij,ij->i',sf[:ni],centers[neighbour]-centers[owner[:ni]])>0)
    closure=np.zeros((nc,3));np.add.at(closure,owner,sf);np.add.at(closure,neighbour,-sf[:ni])
    volume=np.zeros(nc);contribution=np.einsum('ij,ij->i',sf,fc)/3
    np.add.at(volume,owner,contribution);np.add.at(volume,neighbour,-contribution[:ni])
    assert abs(closure).max()<1e-14 and volume.min()>0
    # 模块二：网格几何/诊断控制分离；没有供foamRun使用的物性或初始场。
    mesh=out/'constant/polyMesh';mesh.mkdir(parents=True);(out/'system').mkdir();(out/'0').mkdir()
    def listing(name,kind,values,count):
        with (mesh/name).open('w',encoding='ascii',newline='\n') as f:
            f.write(header(name,kind)+str(count)+'\n(\n')
            for line in values:f.write(line+'\n')
            f.write(')\n// ************************************************************************* //\n')
    listing('points','vectorField',('('+' '.join('%.17g'%a for a in p)+')' for p in P),len(P))
    listing('faces','faceList',(str(n)+'('+' '.join(map(str,f[:n]))+')' for f,n in zip(F,sizes)),len(F))
    listing('owner','labelList',map(str,owner),len(owner));listing('neighbour','labelList',map(str,neighbour),len(neighbour))
    patches=[];offset=ni
    for name,number,kind in [('wall',7838,'wall'),('layerTop',7838,'patch'),('cut',1472,'patch')]:
        patches.append(dict(name=name,nFaces=number,startFace=offset,type=kind));offset+=number
    listing('boundary','polyBoundaryMesh',(p['name']+' { type '+p['type']+'; nFaces '+str(p['nFaces'])+'; startFace '+str(p['startFace'])+'; }' for p in patches),3)
    (out/'system/controlDict').write_text(header('controlDict')+'application checkMesh; startFrom startTime; startTime 0; stopAt endTime; endTime 0; deltaT 1; writeControl timeStep; writeInterval 1;\n')
    (out/'system/fvSchemes').write_text(header('fvSchemes')+'ddtSchemes { default Euler; } gradSchemes { default Gauss linear; } divSchemes { default none; } laplacianSchemes { default Gauss linear corrected; } interpolationSchemes { default linear; } snGradSchemes { default corrected; }\n')
    (out/'system/fvSolution').write_text(header('fvSolution')+'solvers {}\n')
    np.savez_compressed(out/'exact_mesh.npz',points=P,cells=C,faces=F,faceSizes=sizes,owner=owner,neighbour=neighbour,
        originalNodeIDs=old_points+1,originalWedgeRows=original_rows,originalWallFaceIDs=np.array(wall_ids))
    report=dict(passed=True,kind='original-layer-submesh-for-native-diagnostic-only',cells=nc,points=len(P),faces=len(F),internalFaces=ni,
        patches=patches,originalWallFaces=7838,wallAreaM2=float(np.linalg.norm(sf[ni:ni+7838],axis=1).sum()),
        original32LayersRetained=True,protected19LayersRetained=True,originalCoordinatesExactlyPreserved=True,
        vectorAreaMaxClosure=float(abs(closure).max()),volumeM3=float(volume.sum()),minimumCellVolumeM3=float(volume.min()),
        inputSHA256={str(source.relative_to(ROOT)):sha(source),str(selection.relative_to(ROOT)):sha(selection)},
        meshSHA256={p.name:sha(p) for p in mesh.iterdir()},elapsedSeconds=time.monotonic()-start,
        exteriorVolumeMeshed=False,CFDRuns=0,physicalQualification=False,newQualifiedHeatAreaM2=0)
    assert abs(report['wallAreaM2']-.7198603834559664)<1e-13
    shutil.copyfile(__file__,out/'executed_entry.py');save(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['inputSHA256','meshSHA256','patches']}))

if __name__=='__main__':main()
