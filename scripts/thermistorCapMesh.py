"""配对控制只重剖热敏珠两端的小块子午域，切口及其外侧沿用带线算例的实际节点。"""
from collections import Counter, defaultdict

import numpy as np

from thermistorBead import A, D


def edges(triangles):
    return Counter(tuple(sorted((int(a), int(b)))) for t in triangles for a, b in zip(t, np.roll(t, -1)))


def cut_chain(tris, removed, interface, xy):
    border = [e for e, count in edges(tris[removed]).items() if count == 1]
    kept = edges(tris[~removed])
    neighbors = defaultdict(list)
    for a, b in border:
        if (a, b) in kept:
            neighbors[a].append(b)
            neighbors[b].append(a)
    ends = [a for a, b in neighbors.items() if len(b) == 1]
    assert len(ends) == 2 and all(len(n) <= 2 for n in neighbors.values())
    start, = [a for a in ends if a in interface]
    path = [start]
    while True:
        candidates = [a for a in neighbors[path[-1]] if len(path) == 1 or a != path[-2]]
        if not candidates:
            break
        path.append(candidates[0])
    assert len(path) == len(neighbors) and xy[path[-1], 1] == 0
    return path


def patch_triangles(points, chain, front):
    """切口与近壁层接口的每条边固定为一个单元，轴线边可自由剖分。"""
    import gmsh
    assert front[-1] == chain[0]
    loop = chain + front[:-1]
    xy = np.array([points[i] for i in loop])
    cross = np.sum(xy[:, 0]*np.roll(xy[:, 1], -1) - xy[:, 1]*np.roll(xy[:, 0], -1))
    assert cross != 0
    axis = {chain[-1], front[0]}
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.add("B-thermistor-polar-cap")
        geo = gmsh.model.geo
        pts = {}
        for j, i in enumerate(loop):
            step = min(np.linalg.norm(xy[j]-xy[j-1]), np.linalg.norm(xy[j]-xy[(j+1) % len(xy)]))
            pts[i] = geo.addPoint(*points[i], 0, max(.012, min(.08, step)))
        lines, fixed = [], []
        for i, j in zip(loop, loop[1:]+loop[:1]):
            line = geo.addLine(pts[i], pts[j])
            lines.append(line)
            if {i, j} != axis:
                fixed.append(line)
        if cross < 0:
            lines = [-i for i in reversed(lines)]
        face = geo.addPlaneSurface([geo.addCurveLoop(lines)])
        geo.synchronize()
        for line in fixed:
            gmsh.model.mesh.setTransfiniteCurve(line, 2)
        gmsh.model.mesh.generate(2)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        native = {int(gmsh.model.mesh.getNodes(0, p)[0][0]): i for i, p in pts.items()}
        for tag, pos in zip(tags, np.array(coords).reshape(-1, 3)):
            tag = int(tag)
            if tag in native:
                np.testing.assert_array_equal(pos[:2], points[native[tag]])
            else:
                native[tag] = len(points)
                points.append(pos[:2].tolist())
        kinds, _, nodes = gmsh.model.mesh.getElements(2, face)
        assert list(kinds) == [2]
        return np.array([native[int(n)] for n in nodes[0]]).reshape(-1, 3)
    finally:
        gmsh.finalize()


def replace_caps(old_front, new_front, background):
    xy, triangles, far, interface, _ = background
    np.testing.assert_array_equal(xy[interface], old_front)
    points = xy.tolist()
    known = {tuple(x): int(i) for x, i in zip(old_front, interface)}
    new_ids = []
    for p in new_front:
        i = known.get(tuple(p))
        if i is None:
            i = len(points)
            points.append(p.tolist())
        new_ids.append(i)
    removed = np.zeros(len(triangles), bool)
    patches = []
    for side in (-1, 1):
        distance = np.linalg.norm(xy - [side*A/D, 0], axis=1)
        pocket = (distance[triangles] < .45).any(axis=1)
        assert not (pocket & removed).any()
        chain = cut_chain(triangles, pocket, interface, xy)
        j = new_ids.index(chain[0])
        front = new_ids[:j+1] if side == -1 else new_ids[j:][::-1]
        patches.append(patch_triangles(points, chain, front))
        removed |= pocket
    tris = np.vstack([triangles[~removed], *patches])
    assert len(triangles[~removed]) > 0
    return np.asarray(points), tris, far, np.array(new_ids), np.unique(tris)
