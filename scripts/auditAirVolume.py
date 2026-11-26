"""从四面体的真实连接恢复完整边界、原壁相邻单元和体积收支。所有输入编号从零开始。"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from auditSU2HybridTopology import polyvolume

FACES = np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])


def face_sign(faces):
    """从最小编号顶点比较循环方向，不依赖三角面首顶点的选择。"""
    first = np.argmin(faces,axis=1)
    rows = np.arange(len(faces))
    return np.sign(faces[rows,(first+1)%3]-faces[rows,(first-1)%3])


def audit(points, tets, wall, far, box):
    assert np.isfinite(points).all() and len(np.unique(tets))==len(points)
    assert len(np.unique(points, axis=0)) == len(points)
    assert len(np.unique(np.sort(tets, axis=1), axis=0)) == len(tets)
    assert tets.min()==0 and tets.max()==len(points)-1
    corner = points[tets]
    volume = np.einsum("ij,ij->i",corner[:,1]-corner[:,0],np.cross(corner[:,2]-corner[:,0],corner[:,3]-corner[:,0]))/6
    assert np.all(volume>0), "体网包含零体积或负体积单元"
    del corner
    faces = tets[:,FACES].reshape(-1,3)
    keys,inv,count = np.unique(np.sort(faces,axis=1), axis=0, return_inverse=True,return_counts=True)
    inv = inv.ravel()
    assert count.max()==2
    directions = np.bincount(inv,weights=face_sign(faces))
    assert np.all(directions[count==2]==0), "内部面未按相反绕序配对"
    sort = np.argsort(inv,kind="stable")
    head = np.r_[0,np.cumsum(count[:-1])]
    shared = head[count==2]
    graph = coo_matrix((np.ones(len(shared),dtype=np.int8),(sort[shared]//4,sort[shared+1]//4)),shape=(len(tets),len(tets)))
    components = connected_components(graph,directed=False,return_labels=False)
    assert components==1, "空气体网不连通"

    # 空气单元的边界法向朝实体内；STL 的原法向朝实体外。
    wanted = np.concatenate([wall[:,[0,2,1]],far])
    wkeys = np.sort(wanted,axis=1)
    word = np.lexsort(wkeys.T[::-1])
    assert np.array_equal(wkeys[word],keys[count==1]), "实际边界与原壁、外盒的并集不一致"
    single = sort[head[count==1]]
    assert np.array_equal(face_sign(wanted[word]),face_sign(faces[single])), "壁面或外盒绕序错误"
    owners = np.empty(len(wanted),dtype=np.int32)
    owners[word] = single//4
    wall_owner = owners[:len(wall)]
    # 检查每个原面的第四顶点确实位于其正法向一侧。
    triangle = points[wall]
    cross = np.cross(triangle[:,1]-triangle[:,0],triangle[:,2]-triangle[:,0])
    areas = np.linalg.norm(cross,axis=1)/2
    center = points[tets[wall_owner]].mean(axis=1)
    height = 4*np.einsum("ij,ij->i",center-triangle[:,0],cross)/(2*areas)
    assert np.all(height>0), "原壁的相邻体单元位于错误一侧"
    body = polyvolume(points,wall)
    box_volume = np.prod(box[1]-box[0])
    expected = box_volume-body
    error = abs(volume.sum()/expected-1)
    assert body>0 and expected>0 and error<=2e-10
    assert abs(polyvolume(points,far)/box_volume-1)<=2e-10
    return wall_owner,areas,dict(nodes=len(points),tetrahedra=len(tets),wallFaces=len(wall),outerFaces=len(far),
        internalFaces=int(np.sum(count==2)),fluidComponents=int(components),wallAreaM2=float(areas.sum()),
        bodyVolumeM3=float(body),boxVolumeM3=float(box_volume),airVolumeM3=float(volume.sum()),relativeVolumeError=float(error),
        minTetVolumeM3=float(volume.min()),minimumWallOppositeHeightM=float(height.min()))
