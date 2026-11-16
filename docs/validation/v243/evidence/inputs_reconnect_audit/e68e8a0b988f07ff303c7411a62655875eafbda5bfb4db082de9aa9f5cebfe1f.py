"""模块职责：将现有四面体外流网格裁为共形平面出口，保留原壁面和母网格来源。"""
from collections import defaultdict
from pathlib import Path
import numpy as np


def triangulate(poly):
    """模块一：共享多边形采用最小全局编号扇形，邻单元得到相同对角线。"""
    poly = list(dict.fromkeys(poly))
    if len(poly) < 3:
        return []
    k = poly.index(min(poly))
    poly = poly[k:] + poly[:k]
    return [[poly[0], poly[j], poly[j+1]] for j in range(1, len(poly)-1)]


def clip_tetrahedra(xyz, cells, markers, cut_x):
    """模块二：保留x>=cut_x的多面体；原三角面线性裁切，部分单元作内部星形剖分。

    不移动任何原节点。分类容差仅识别共面舍入差，不能用来修补不在平面的壁环。
    返回的稀疏线性权重只用于几何/已存研究场溯源，不自动授予CFD重启动资格。
    """
    xyz = np.asarray(xyz, dtype=float)
    cells = np.asarray(cells, dtype=np.int64)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or cells.ndim != 2 or cells.shape[1] != 4:
        raise ValueError('CLIP_REQUIRES_TETRA_MESH')
    if not np.isfinite(xyz).all() or not np.isfinite(cut_x):
        raise ValueError('CLIP_NONFINITE_COORDINATE')
    if np.any(cells < 0) or np.any(cells >= len(xyz)):
        raise ValueError('CLIP_NODE_ID')
    tol = 64*np.finfo(float).eps*max(1., float(np.max(abs(xyz))), abs(cut_x))
    signed = xyz[:, 0]-cut_x
    side = np.where(abs(signed) <= tol, 0, np.sign(signed)).astype(int)
    coords = list(xyz.copy())
    origins = [{int(i): 1.} for i in range(len(xyz))]
    edges = {}

    def intersection(a, b):
        if side[a] == 0: return int(a)
        if side[b] == 0: return int(b)
        key = tuple(sorted((int(a), int(b))))
        if key not in edges:
            i, j = key
            t = (cut_x-xyz[i, 0])/(xyz[j, 0]-xyz[i, 0])
            if not 0 < t < 1: raise ValueError('CLIP_EDGE_NOT_STRADDLING')
            q = xyz[i] + t*(xyz[j]-xyz[i]); q[0] = cut_x
            edges[key] = len(coords); coords.append(q)
            origins.append({i: 1-t, j: t})
        return edges[key]

    def polygon(face):
        result = []
        for a, b in zip(face, np.roll(face, -1)):
            if side[a] >= 0: result.append(int(a))
            if side[a]*side[b] < 0: result.append(intersection(a, b))
        return list(dict.fromkeys(result))

    new_cells, parent, outlet = [], [], []
    full = np.all(side[cells] >= 0, axis=1)
    active = np.any(side[cells] > 0, axis=1)
    full &= active
    for i in np.flatnonzero(full):
        new_cells.append(cells[i].tolist()); parent.append(int(i))
    partial = np.flatnonzero(active & ~full)
    face_indices = ((0, 1, 2), (0, 3, 1), (1, 3, 2), (2, 3, 0))
    for i in partial:
        tet = cells[i]
        polygons = [polygon(tet[list(f)]) for f in face_indices]
        vertices = sorted(set(v for p in polygons for v in p))
        cap = [v for v in vertices if abs(coords[v][0]-cut_x) <= tol]
        if len(cap) < 3: raise ValueError('CLIP_MISSING_CAP')
        center = np.mean([coords[v] for v in cap], axis=0)
        cap.sort(key=lambda v: np.arctan2(coords[v][2]-center[2], coords[v][1]-center[1]))
        cap_triangles = triangulate(cap); outlet.extend(cap_triangles)
        polygons.append(cap)
        if len(vertices) == 4:
            new_cells.append(vertices); parent.append(int(i))
        else:
            c = len(coords); coords.append(np.mean([coords[v] for v in vertices], axis=0))
            weights = defaultdict(float)
            for v in vertices:
                for p, w in origins[v].items(): weights[p] += w/len(vertices)
            origins.append(dict(weights))
            for poly in polygons:
                for tri in triangulate(poly):
                    new_cells.append([c]+tri); parent.append(int(i))
    # 整单元有面恰在切平面时，该面也是新出口；内部共享共面面随后由拓扑审核拒绝。
    for i in np.flatnonzero(full):
        for f in face_indices:
            face = cells[i, list(f)]
            if np.all(side[face] == 0): outlet.append(face.tolist())
    new_markers = {}
    for tag, faces in markers.items():
        if tag == 'outlet': raise ValueError('CLIP_OUTLET_TAG_COLLISION')
        result = []
        for face in faces:
            if np.all(side[face] == 0): continue
            result.extend(triangulate(polygon(face)))
        new_markers[tag] = np.asarray(result, dtype=np.int64).reshape(-1, 3)
    new_markers['outlet'] = np.asarray(outlet, dtype=np.int64).reshape(-1, 3)
    if not len(new_cells) or not len(outlet): raise ValueError('CLIP_EMPTY_REGION_OR_OUTLET')
    coords = np.asarray(coords); new_cells = np.asarray(new_cells, dtype=np.int64)
    used = np.unique(new_cells); remap = np.full(len(coords), -1, dtype=np.int64); remap[used] = np.arange(len(used))
    result_xyz = coords[used]; new_cells = remap[new_cells]
    new_markers = {k: remap[v] for k, v in new_markers.items()}
    if any(np.any(v < 0) for v in new_markers.values()): raise ValueError('CLIP_UNUSED_BOUNDARY_NODE')
    delta = result_xyz[new_cells[:, 1:]]-result_xyz[new_cells[:, 0, None]]
    det = np.einsum('ij,ij->i', delta[:, 0], np.cross(delta[:, 1], delta[:, 2]))
    if np.any(det == 0): raise ValueError('CLIP_ZERO_VOLUME')
    neg = det < 0
    new_cells[neg, 1], new_cells[neg, 2] = new_cells[neg, 2].copy(), new_cells[neg, 1].copy()
    source_ids = np.zeros((len(used), 4), dtype=np.int64); weights = np.zeros((len(used), 4))
    for i, old in enumerate(used):
        if len(origins[old]) > 4: raise ValueError('CLIP_NONLOCAL_ORIGIN')
        for j, (p, w) in enumerate(sorted(origins[old].items())):
            source_ids[i, j] = p; weights[i, j] = w
    old_ids = np.where(used < len(xyz), used, -1)
    return dict(xyz=result_xyz, cells=new_cells, markers=new_markers, motherCell=np.asarray(parent),
                originalNode=old_ids, sourceNode=source_ids, sourceWeight=weights,
                cutX=float(cut_x), classificationTolerance=tol, partialMotherCells=partial,
                fullMotherCells=np.flatnonzero(full))


def interpolate_saved_field(mesh, values):
    """模块三：研究场按同一原单元仿射坐标插值，原节点直接复制以保留原始位形。"""
    values = np.asarray(values)
    result = np.sum(values[mesh['sourceNode']]*mesh['sourceWeight'][:, :, None], axis=1)
    old = mesh['originalNode'] >= 0
    result[old] = values[mesh['originalNode'][old]]
    return result


def write_su2(path, mesh):
    """模块四：写真实可运行SU2网格；保留wall/farfield/outlet的独立边界身份。"""
    with Path(path).open('w', encoding='ascii', newline='\n') as f:
        f.write('NDIME= 3\nNELEM= %d\n' % len(mesh['cells']))
        for i, row in enumerate(mesh['cells']): f.write('10 %d %d %d %d %d\n' % tuple(list(row)+[i]))
        f.write('NPOIN= %d\n' % len(mesh['xyz']))
        for i, row in enumerate(mesh['xyz']): f.write('%.17g %.17g %.17g %d\n' % tuple(list(row)+[i]))
        f.write('NMARK= %d\n' % len(mesh['markers']))
        for tag, faces in mesh['markers'].items():
            f.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n' % (tag, len(faces)))
            for row in faces: f.write('5 %d %d %d\n' % tuple(row))
