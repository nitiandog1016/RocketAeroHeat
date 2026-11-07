# -*- coding: utf-8 -*-
"""原生Euler壁邻采样位置合同；记录实际体单元位置，禁止冒充边界层外缘。"""
import numpy as np
from readEilmerNative import need,sha
from hybridExternalState import fingerprint


def map_original_wall(surface, grid, wall_mapping, cell_positions):
    """模块一：按原三角/四边面连接双射映射，保留原面顺序及唯一邻单元。"""
    vertices=np.asarray(surface['vertices'],float);triangles=np.asarray(surface['facesOneBased'],int)-1
    need(vertices.ndim==2 and vertices.shape[1]==3 and len(vertices)<=len(grid['points']),'原壁顶点尺寸错误')
    need(np.array_equal(vertices,grid['points'][:len(vertices)]),'原壁顶点身份改变')
    need(triangles.ndim==2 and triangles.shape[1]==3 and triangles.min()>=0 and triangles.max()<len(vertices),'原壁面连接错误')
    need(surface['faceIDs']==list(range(1,len(triangles)+1)),'原面ID不连续或重排')
    b=grid['boundaries']['wall'];source=wall_mapping['originalNodes'];rows=np.asarray(wall_mapping['nativeRows'],int)
    need(len(rows)==len(source)==len(b['ids']) and set(rows)==set(range(len(rows))),'原生壁面映射非双射')
    table={};native_triangles={}
    for nodes,row in zip(source,rows):
        fid=int(b['ids'][row]);native=grid['faces'][fid];native=native[native>=0]
        need(set(nodes)==set(native),'原生面与源连接不相容')
        need(len(nodes) in (3,4),'不支持的壁面类型')
        pieces=[nodes] if len(nodes)==3 else [[nodes[0],nodes[1],nodes[2]],[nodes[0],nodes[2],nodes[3]]]
        for piece in pieces:
            key=tuple(sorted(piece));need(key not in table,'重复原壁三角面')
            table[key]=(int(row),fid,int(b['owners'][row]));native_triangles[key]=native
    need(len(table)==len(triangles),'遗漏原壁三角面')
    keys=[tuple(sorted(t)) for t in triangles]
    need(len(set(keys))==len(keys) and set(keys)==set(table),'原面缺失、重复或三角剖分改变')
    positions=np.asarray(cell_positions,float)
    need(positions.shape==(len(grid['cells']),3) and np.isfinite(positions).all(),'原生采样单元位置错误')
    ids=np.array([table[k] for k in keys],int);owners=ids[:,2];centers=vertices[triangles].mean(1)
    # 原生boundary signs确定流体域外向。原箭体法向与控制面文件方向可不同，不能猜号。
    normals=[]
    for key,(row,fid,owner) in zip(keys,ids):
        xyz=grid['points'][native_triangles[key]];vec=np.zeros(3)
        for j in range(1,len(xyz)-1):vec+=np.cross(xyz[j]-xyz[0],xyz[j+1]-xyz[0])/2
        need(np.linalg.norm(vec)>0,'原生退化壁面')
        normals.append(-int(b['signs'][row])*vec/np.linalg.norm(vec))
    normals=np.asarray(normals);delta=positions[owners]-centers
    distance=np.einsum('ij,ij->i',delta,normals)
    need(np.all(distance>0),'壁邻采样位置未在流体侧')
    return dict(originalFaceIDs=surface['faceIDs'],nativeFaceIDsZeroBased=ids[:,1].tolist(),
        nativeOwnerCellsZeroBased=owners.tolist(),wallCentersM=centers.tolist(),
        samplePositionsM=positions[owners].tolist(),normalDistanceM=distance.tolist(),
        tangentialOffsetM=np.linalg.norm(delta-distance[:,None]*normals,axis=1).tolist(),
        sampleKind='native-finite-volume-cell-average-at-wall-neighbour',
        positionKind='native-exported-cell-position',mappingMethod='exact-original-connectivity-bijection')


def candidate(frame,surface,grid,mapping,raw_fields):
    """模块二：状态值逐项核对所指单元；采样位置完善不增加物理/外缘资格。"""
    need(frame.get('schema')=='RAH.HybridExternalState.v1' and frame.get('algebraicInterfacePassed'),'缺少同EOS已审核原生交换帧')
    need(frame.get('surfaceSHA256')==fingerprint(surface),'原表面SHA不相容')
    need(bool(frame.get('sourceFilesSHA256')),'缺少实际源文件')
    for path,digest in frame['sourceFilesSHA256'].items():need(sha(path)==digest,'源帧文件已变')
    need(frame.get('faceIDs')==surface['faceIDs'],'帧与原面身份不符')
    positions=np.column_stack([raw_fields['pos.'+axis] for axis in 'xyz'])
    record=map_original_wall(surface,grid,mapping['wall'],positions)
    owners=np.asarray(record['nativeOwnerCellsZeroBased']);fields=frame['fields']
    for key,raw in [('pressurePa','p'),('temperatureK','T'),('densityKgm3','rho'),('internalEnergyJkg','e')]:
        actual=np.asarray(fields[key]);expected=raw_fields[raw][owners]
        need(actual.shape==expected.shape and np.all(np.isfinite(actual)) and
             np.all(abs(actual-expected)<=2e-10*np.maximum(1,abs(expected))),'所指单元状态不一致: '+key)
    velocity=np.column_stack([raw_fields['vel.'+axis][owners] for axis in 'xyz'])
    need(np.array_equal(velocity,np.asarray(fields['velocityMps'])),'采样三维速度与原帧不一致')
    return dict(schema='RAH.HybridEdgeCandidate.v1',sourceFrameSchema=frame['schema'],
        surfaceSHA256=frame['surfaceSHA256'],conditionSHA256=frame['conditionSHA256'],modelSHA256=frame['modelSHA256'],
        sourceFilesSHA256=frame['sourceFilesSHA256'],snapshot=frame['snapshot'],
        originalAreaM2=frame['originalAreaM2'],faceIDs=frame['faceIDs'],fields=fields,sampling=record,
        temperatureIdentity='Euler wall-neighbour fluid temperature; not Te/Tr/solid Tw',
        algebraicInterfacePassed=True,samplingIdentityPassed=True,attachedEdgeQualified=False,
        engineeringHeatReady=False,physicalQualification=False,
        pendingEvidence=['spatial and farfield sensitivity of local p/T/U/h0/s',
            'streamline and shock/entropy-layer provenance','attached-boundary-layer applicability',
            'same-model energy/entropy/transport or validated cross-model coupling','matched public thermal validation'],
        wallHeatFluxWm2=[None]*len(owners),wallTemperatureSolvedK=[None]*len(owners))
