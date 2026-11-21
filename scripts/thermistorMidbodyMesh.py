"""只细化光滑热敏珠中段的子午三角形，原壁、棱柱层接口与远区边保持原位。"""
import numpy as np

from thermistorBead import A, D


def wall_distance(points, wall):
    """量取点至实际分段直线母线的最短距离，所有坐标均已除以横径。"""
    distance = np.full(len(points), np.inf)
    for start, end in zip(wall, wall[1:]):
        direction = end - start
        offset = points - start
        fraction = np.clip((offset @ direction) / (direction @ direction), 0, 1)
        distance = np.minimum(distance, np.linalg.norm(offset - fraction[:, None]*direction, axis=1))
    return distance


def refine_midbody(background, wall):
    """共享边只建一个中点。两次局部二分允许中段变细，邻接三角形共用同一切分。"""
    xy, triangles, far, interface, _ = background
    for _ in range(2):
        pairs = triangles[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2)
        edges, inverse, counts = np.unique(np.sort(pairs, axis=1), axis=0,
                                           return_inverse=True, return_counts=True)
        ends = xy[edges]
        centers = ends.mean(axis=1)
        distance = wall_distance(centers, wall)
        axial = abs(centers[:, 0]) / (A/D)
        strength = np.clip((.85-axial)/.25, 0, 1) * np.clip((.16-distance)/.06, 0, 1)
        split = (counts == 2) & (np.linalg.norm(ends[:, 1]-ends[:, 0], axis=1)*strength > .015)
        mid = np.full(len(edges), -1, dtype=int)
        mid[split] = len(xy) + np.arange(np.count_nonzero(split))
        xy = np.vstack((xy, centers[split]))
        children = []
        for vertices, mids in zip(triangles, mid[inverse].reshape(-1, 3)):
            n = np.count_nonzero(mids >= 0)
            if n == 0:
                children.append(vertices)
            elif n == 3:
                a, b, c = vertices
                u, v, w = mids
                children.extend(((a, u, w), (u, b, v), (w, v, c), (u, v, w)))
            else:
                # 旋转顶点编号，让第一条边成为唯一已分边，或让最后一条成为唯一未分边。
                i = int(np.flatnonzero(mids >= 0)[0]) if n == 1 else (int(np.flatnonzero(mids < 0)[0])+1) % 3
                a, b, c = np.roll(vertices, -i)
                u, v, _ = np.roll(mids, -i)
                if n == 1:
                    children.extend(((a, u, c), (u, b, c)))
                else:
                    children.append((u, b, v))
                    if np.linalg.norm(xy[u]-xy[c]) <= np.linalg.norm(xy[a]-xy[v]):
                        children.extend(((a, u, c), (u, v, c)))
                    else:
                        children.extend(((a, u, v), (a, v, c)))
        triangles = np.asarray(children)
    return xy, triangles, far, interface, np.unique(triangles)
