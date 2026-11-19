"""核对导入的 CM-011 原面，用实际环坐标生成固定低速试件的求解网格。"""
import argparse
import json
from pathlib import Path
import shutil
import uuid

import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree

from cm011Inclined import prepare
from cm011Tunnel import D, ROOT, meridian, surface, write_json
from runEmbeddedSU2 import sha

TOL_M = 5e-10
FILES = ("body.json","mesh.su2","case.cfg","nodes.npz","geometry.json","initial_flow.csv")


def oriented_keys(faces):
    # 将最小顶点移到第一列；循环换序保持方向，反向面会产生不同的键。
    order = (np.argmin(faces,axis=1)[:,None]+np.arange(3))%3
    return np.take_along_axis(faces,order,axis=1)


def match_surface(vertices,faces):
    line,bands = meridian(2)
    ref = surface(line,bands,2)
    expected = np.array(ref["V"])
    ref_faces = np.array(ref["F"])-1
    if (vertices.shape!=expected.shape or faces.shape!=ref_faces.shape or not np.isfinite(vertices).all() or
            not np.array_equal(faces,faces.astype(int)) or faces.min()<1 or faces.max()>len(vertices)):
        raise ValueError("导入表面的节点或三角面与 CM-011 缩比试件不匹配")
    idx = faces.astype(int)-1
    dist,vertex_id = cKDTree(expected).query(vertices)
    if max(dist)>TOL_M or len(np.unique(vertex_id))!=len(expected):
        raise ValueError("CM-011 节点匹配失败，请核对尺度、朝向和局部几何")
    lookup = {tuple(key):i for i,key in enumerate(oriented_keys(ref_faces))}
    face_id = np.array([lookup.get(tuple(key),-1) for key in oriented_keys(vertex_id[idx])])
    if np.any(face_id<0) or len(np.unique(face_id))!=len(ref_faces):
        raise ValueError("导入面与参考三角面无法同向一一对应")
    rows = np.array(ref["nodeRow"])[vertex_id]
    radii = np.hypot(vertices[:,1],vertices[:,2])
    order = np.argsort(vertex_id)
    actual = np.column_stack([np.bincount(rows[order],weights=v[order])/np.bincount(rows)
                              for v in (vertices[:,0],radii)])/D
    if np.any(actual[[0,-1],1]>TOL_M/D):
        raise ValueError("两个极点偏离试件轴线")
    actual[[0,-1],1] = 0.
    ring_error = np.hypot(vertices[:,0]-D*actual[rows,0],radii-D*actual[rows,1])
    body = dict(V=vertices.tolist(),F=(idx+1).tolist(),nodeRow=rows.tolist(),
                regionCode=np.array(ref["regionCode"])[face_id].tolist(),regionNames=ref["regionNames"],
                meridian=actual.tolist(),diameterM=D,grid="fine",extentD=12.)
    reg = dict(vertexToleranceM=TOL_M,maxVertexDistanceM=float(max(dist)),maxRingDeviationM=float(max(ring_error)),
               referenceVertexOneBased=(vertex_id+1).tolist(),referenceFaceOneBased=(face_id+1).tolist())
    return body,reg


def prepare_import(folder,mat,source,source_hash,scale):
    if not np.isfinite(scale) or scale<=0 or sha(source)!=source_hash:
        raise ValueError("STL 尺度无效或文件自导入后发生改变")
    arrays=loadmat(mat,variable_names=("V","F"))
    body,reg=match_surface(arrays["V"],arrays["F"])
    reg.update(sourceFile=source.name,sourceSHA256=source_hash,scaleToM=scale,inputArraySHA256=sha(mat))
    prepare(folder,body,reg)
    write_json(folder/"prepared.json",{n:sha(folder/n) for n in FILES})


def instantiate(template,folder,expected):
    if sha(template/"prepared.json")!=expected:
        raise ValueError("体网准备回执与当前确认的网格不匹配")
    manifest=json.loads((template/"prepared.json").read_bytes())
    if set(manifest)!=set(FILES):
        raise ValueError("体网准备回执缺少输入")
    folder.mkdir(parents=True,exist_ok=False)
    for n in FILES:
        shutil.copyfile(template/n,folder/n)
        if sha(folder/n)!=manifest[n]:
            raise ValueError("输入发生变化："+n)
    shutil.copyfile(template/"prepared.json",folder/"prepared.json")
    if sha(folder/"prepared.json")!=expected:
        raise ValueError("体网准备回执已更改")
    sources=[Path(__file__),Path(__file__).with_name("cm011Inclined.py"),folder/"prepared.json",folder/"body.json",
             folder/"nodes.npz",folder/"geometry.json"]
    write_json(folder/"engine_request.json",dict(workDirectory=str(folder),configFile="case.cfg",meshFile="mesh.su2",
               maximumSeconds=720,threads=1,jobToken=uuid.uuid4().hex,
               sourceFiles=[p.resolve().relative_to(ROOT).as_posix() for p in sources]))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="action",required=True)
    prep=sub.add_parser("prepare")
    for name in ("directory","mat","source"): prep.add_argument(name,type=Path)
    prep.add_argument("sha")
    prep.add_argument("scale",type=float)
    clone=sub.add_parser("instantiate")
    clone.add_argument("template",type=Path)
    clone.add_argument("directory",type=Path)
    clone.add_argument("sha")
    args=parser.parse_args()
    folder=args.directory.resolve()
    folder.relative_to(ROOT/"outputs")
    if args.action=="prepare": prepare_import(folder,args.mat,args.source,args.sha,args.scale)
    else: instantiate(args.template,folder,args.sha)
