"""用已有闭合表面作为固定内边界，生成包围盒内的四面体空气域。坐标为米。"""
import sys

import numpy as np

import installGmshRuntime


def mesh_air(vertices, wall):
    installGmshRuntime.main()
    runtime = installGmshRuntime.ROOT / "tools/runtime/gmsh" / installGmshRuntime.VERSION
    sys.path.insert(0, str(runtime))
    import gmsh

    length = np.ptp(vertices, axis=0).max()
    lo, hi = vertices.min(axis=0)-length/2, vertices.max(axis=0)+length/2
    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.model.add("closed-surface-air")
        gmsh.option.setNumber("General.NumThreads", 1)
        # 先建立盒壁，再添加离散内壳，两种几何使用不同的实体编号。
        box = gmsh.model.occ.addBox(*lo, *(hi-lo))
        gmsh.model.occ.synchronize()
        sides = [tag for _,tag in gmsh.model.getBoundary([(3,box)], oriented=False)]
        gmsh.model.occ.remove([(3,box)], recursive=False)
        gmsh.model.occ.synchronize()
        body = gmsh.model.addDiscreteEntity(2, max(sides)+1)
        ids = np.arange(1,len(vertices)+1)
        gmsh.model.mesh.addNodes(2, body, ids, vertices.ravel())
        gmsh.model.mesh.addElementsByType(body, 2, np.arange(1,len(wall)+1), (wall+1).ravel())
        boundary = gmsh.model.geo.addSurfaceLoop(sides)
        hole = gmsh.model.geo.addSurfaceLoop([body])
        air = gmsh.model.geo.addVolume([boundary,hole])
        gmsh.model.geo.synchronize()
        for name,value in (("Mesh.MeshOnlyEmpty",1), ("Mesh.Renumber",0), ("Mesh.Algorithm3D",1),
                           ("Mesh.MeshSizeExtendFromBoundary",0), ("Mesh.MeshSizeMax",length/5)):
            gmsh.option.setNumber(name,value)
        gmsh.model.mesh.generate(3)

        tags, xyz, _ = gmsh.model.mesh.getNodes()
        order = np.argsort(tags)
        tags = tags[order]
        points = xyz.reshape(-1,3)[order]
        if not np.array_equal(tags[:len(ids)],ids) or not np.array_equal(points[:len(ids)],vertices):
            raise ValueError("建体后原壁节点移动或编号改变。")
        def elements(dim, entity, kind, width):
            kinds, numbers, nodes = gmsh.model.mesh.getElements(dim,entity)
            if list(kinds)!=[kind]:
                raise ValueError(f"实体 {entity} 出现未登记的单元类型 {kinds}。")
            raw = nodes[0].reshape(-1,width)
            idx = np.searchsorted(tags,raw)
            assert np.array_equal(tags[idx],raw)
            return idx.astype(np.int32),numbers[0]

        tets, native_ids = elements(3,air,4,4)
        if len(tets)>1200000 or not len(tets):
            raise ValueError("四面体数量超出本次建体预算。")
        far = np.concatenate([elements(2,s,2,3)[0] for s in sides])
        outer_points = points[far]
        normal = np.cross(outer_points[:,1]-outer_points[:,0],outer_points[:,2]-outer_points[:,0])
        flip = np.einsum("ij,ij->i",normal,outer_points.mean(axis=1)-(lo+hi)/2)<0
        far[flip] = far[flip][:,[0,2,1]]
        quality = gmsh.model.mesh.getElementQualities(native_ids,"minSICN")
        assert len(quality) == len(tets) and np.all(np.isfinite(quality) & (quality > 0))
        return points, tets, far, np.array([lo, hi]), quality
    finally:
        gmsh.finalize()
