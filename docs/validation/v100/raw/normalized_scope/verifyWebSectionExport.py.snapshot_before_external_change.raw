"""独立核对浏览器 CSV：来源、三角形交点、轮廓拓扑/角度及圆柱周长。"""
import argparse
import csv
import json
import math
from pathlib import Path


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def inside_triangle(point, triangle):
    a, b, c = triangle
    u, v, w = sub(b, a), sub(c, a), sub(point, a)
    uu, vv, uv, wu, wv = dot(u, u), dot(v, v), dot(u, v), dot(w, u), dot(w, v)
    denominator = uu * vv - uv * uv
    beta = (vv * wu - uv * wv) / denominator
    gamma = (uu * wv - uv * wu) / denominator
    residual = math.sqrt(sum((w[j] - beta * u[j] - gamma * v[j]) ** 2 for j in range(3)))
    assert min(beta, gamma, 1 - beta - gamma) >= -1e-7
    assert residual < 1e-8


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file", type=Path)
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--field", default="pressure")
    parser.add_argument("--radius", type=float)
    args = parser.parse_args()
    data = json.loads(args.json_file.read_text(encoding="utf-8-sig"))
    field = next(f for f in data["fields"] if f["key"] == args.field)
    assert not field.get("history"), "本校核只针对单时刻数据，不把末站当任意历史站。"
    with args.csv_file.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    assert "角度状态" in rows[0], "必须核对当前带角度状态的 CSV，不能复用旧版导出。"
    extent = max(max(v[j] for v in data["vertices"]) - min(v[j] for v in data["vertices"]) for j in range(3))
    tolerance = max(extent * 1e-10, max(abs(q) for v in data["vertices"] for q in v) * 128 * 2.220446049250313e-16)
    nodes, node_groups, groups, seen_faces = [], [], {}, set()
    expected_x, expected_trace = rows[0]["截面 X / m"], rows[0]["网格环重合时取值侧"]
    assert expected_trace in ("+X", "-X")
    value_key = f'{field["label"]} / {field["unit"]}'
    perimeter, given, missing = 0.0, 0, 0
    for row in rows:
        index = int(row["面元编号"]) - 1
        assert index not in seen_faces, "同一三角面不能重复导出同一平面的非退化交线。"
        seen_faces.add(index)
        assert row["截面 X / m"] == expected_x and row["网格环重合时取值侧"] == expected_trace
        expected, state = field["values"][index], field["status"][index]
        if state == 0 or expected is None:
            assert row[value_key] == "" and row["状态"] == "未求解"
            missing += 1
        else:
            assert math.isclose(float(row[value_key]), expected, rel_tol=1e-14, abs_tol=1e-14)
            assert row["状态"] == ["未求解", "工程计算", "给定条件"][state]
            given += state == 2
        assert row["算法"] == field["methodNames"][field["methodIndex"][index]]
        expected_scope = field["scopeNames"][field["scopeIndex"][index]] if field.get("scopeNames") else field["scope"]
        assert row["适用范围"] == expected_scope
        x = float(row["截面 X / m"])
        a = [x, float(row["A_Y / m"]), float(row["A_Z / m"])]
        b = [x, float(row["B_Y / m"]), float(row["B_Z / m"])]
        triangle = [data["vertices"][i] for i in data["faces"][index]]
        inside_triangle(a, triangle)
        inside_triangle(b, triangle)
        delta = sub(a, b)
        length = math.sqrt(dot(delta, delta))
        assert math.isclose(length, float(row["交线长度 / m"]), rel_tol=1e-12, abs_tol=1e-14)
        perimeter += length
        group = int(row["轮廓编号"])
        endpoint_ids = []
        # 逐点实际距离合并，不复用 JS 的拓扑键或空间桶；共享点不得被拆到不同轮廓。
        for point in (a, b):
            matches = [j for j, other in enumerate(nodes) if math.sqrt(dot(sub(point, other), sub(point, other))) <= tolerance]
            if matches:
                node = matches[0]
                assert node_groups[node] == group, "共享交点被错误拆为多个轮廓。"
            else:
                node = len(nodes)
                nodes.append(point)
                node_groups.append(group)
            endpoint_ids.append(node)
        groups.setdefault(group, []).append((row, a, b, length, endpoint_ids))
    for entries in groups.values():
        weight = sum(e[3] for e in entries)
        center = [float(expected_x)] + [sum((e[1][j] + e[2][j]) * e[3] / 2 for e in entries) / weight for j in (1, 2)]
        degree = {}
        for entry in entries:
            for node in entry[4]:
                degree[node] = degree.get(node, 0) + 1
        closed, branched = all(d == 2 for d in degree.values()), any(d > 2 for d in degree.values())
        if args.radius is not None:
            assert closed and not branched, "圆柱截面应为闭合无分支轮廓。"
        for row, a, b, length, _ in entries:
            av, bv = sub(a, center), sub(b, center)
            delta = sub(b, a)
            t = min(1, max(0, -dot(av, delta) / length ** 2))
            center_distance = math.sqrt(sum((av[j] + t * delta[j]) ** 2 for j in (1, 2)))
            radial = abs(av[1] * bv[2] - av[2] * bv[1]) <= tolerance * (math.hypot(*av[1:]) + math.hypot(*bv[1:]))
            reason = "轮廓存在分支" if branched else "轮廓开放" if not closed else "交线经过角度中心" if center_distance <= tolerance else "交线为径向线，角度跨度退化" if radial else "有效周向角"
            assert row["角度状态"] == reason
            if reason != "有效周向角":
                assert row["角度A / deg"] == "" and row["角度B / deg"] == ""
            else:
                for endpoint, key in ((a, "角度A / deg"), (b, "角度B / deg")):
                    expected_angle = math.degrees(math.atan2(endpoint[2] - center[2], endpoint[1] - center[1])) % 360
                    difference = abs(float(row[key]) - expected_angle) % 360
                    assert min(difference, 360 - difference) < 1e-9
    if args.radius is not None:
        assert math.isclose(perimeter, 2 * math.pi * args.radius, rel_tol=0.005)
    print(f"Actual section CSV passed: {len(rows)} segments, {given} prescribed, {missing} missing; {len(groups)} contours; perimeter {perimeter:.12g} m; original values/methods/scopes, triangle intersections, component topology and angles verified.")


if __name__ == "__main__":
    main()
