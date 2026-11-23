"""用十进制、高精度投影和截面解析积分独立复核，不调用候选实现。"""
import argparse
import csv
from decimal import Decimal as D, getcontext
import hashlib
import io
import json
from pathlib import Path
import time
from zipfile import ZipFile

getcontext().prec = 50
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EVIDENCE = HERE/"evidence"
FIELDS = {"constant": (0, 0, 0), "x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1), "oblique": (1, -2, 3)}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def projected(displacement, area, g, k):
    # 拆成平行及垂直连线两个有向面积向量，算受限结果时从两点通量的区间投影原通量。
    beta = D(dot(displacement, area))/D(dot(displacement, displacement))
    normal_part = [a-beta*d for a, d in zip(area, displacement)]
    base = -k*beta*dot(g, displacement)
    original = base-k*dot(g, normal_part)
    bounded = max(base-abs(base), min(base+abs(base), original))
    return base, bounded, k*beta


def equal(value, expected, scale):
    value = D(str(value))
    assert abs(value-expected) <= D("3e-12")*scale, (value, expected, scale)


def verify(result_file, raw):
    started = time.perf_counter()
    result = json.loads(result_file.read_bytes())
    receipt = json.loads((raw/"receipt.json").read_bytes())
    assert receipt["program"] == sha((HERE/"check_diffusion.py").read_bytes())
    for name, digest in receipt["outputs"].items():
        path = result_file if name == "check.json" else raw/name
        assert sha(path.read_bytes()) == digest
    archive = ROOT/"docs/validation/v326/evidence/receipts.zip"
    assert sha(archive.read_bytes()) == "37b679a7181e5d1a45855242f7caef298f8387a1a47e4ebd6fa882562ca9dca5"
    with ZipFile(archive) as z:
        nodes = {int(r["id"]): tuple(D(r[k]) for k in ("xD", "yD", "zD"))
                 for r in csv.DictReader(io.StringIO(z.read("first_nodes.csv").decode()))}
        original = csv.reader(io.StringIO(z.read("outputs/v326_validation/detail01/scalar_edges.csv").decode()))
        native = {int(r[0]): (int(r[1]), int(r[2]), list(map(D, r[3:6]))) for r in original}
        spec = json.loads(z.read("outputs/v325_validation/prep02/prepared/prepared.json"), parse_float=D)
    k = 1/spec["condition"]["pe"]
    with (raw/"star.csv").open(encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    expected = {(field, edge) for field in FIELDS for edge in native}
    assert len(rows) == len(expected) and {(r["field"], int(r["edge"])) for r in rows} == expected
    actual_ledger = {name: D(0) for name in FIELDS}
    max_error = D(0)
    for r in rows:
        i, j, normal = native[int(r["edge"])]
        g = FIELDS[r["field"]]
        dx = [b-a for a, b in zip(nodes[i], nodes[j])]
        base, bounded, diag = projected(dx, normal, g, k)
        exact = -k*dot(normal, g)
        scale = k*dot(normal, normal).sqrt()*(D(dot(g, g)).sqrt() or 1)
        assert i == int(r["i"]) and j == int(r["j"])
        for name, calc in (("twoPoint", base), ("limited", bounded), ("reverse", -bounded),
                           ("original", exact), ("exact", exact), ("scale", scale)):
            equal(r[name], calc, scale)
        equal(r["diagonal"], diag, diag)
        error = abs(bounded-exact)/scale
        equal(r["normalizedError"], error, D(1))
        if error > max_error:
            max_error, worst_edge, worst_field = error, int(r["edge"]), r["field"]
        sign = 1 if i == result["star"]["pointID"] else -1
        assert int(r["targetSign"]) == sign
        actual_ledger[r["field"]] += sign*(bounded-exact)
    assert (worst_edge, worst_field) == (result["star"]["worst"]["edge"], result["star"]["worst"]["field"])
    with (raw/"cuts.csv").open(encoding="ascii") as stream:
        cuts = list(csv.DictReader(stream))
    combinations = {(s, a, n, f) for s in (D(0), D(".5")) for a in range(3) for n in (3, 9, 27) for f in FIELDS}
    assert {(D(r["shear"]), int(r["axis"]), int(r["divisions"]), r["field"]) for r in cuts} == combinations
    assert len(cuts) == len(combinations) == result["cutRows"]
    for r in cuts:
        shear, axis, field = D(r["shear"]), int(r["axis"]), r["field"]
        g = FIELDS[field]
        line = [(1,0,0), (shear,1,0), (0,0,1)][axis]
        whole = [(1,-shear,0), (0,1,0), (0,0,1)][axis]
        base, limited, _ = projected(line, whole, g, k)
        full = -k*dot(whole, g)
        scale = k*D(dot(whole,whole)).sqrt()*(D(dot(g,g)).sqrt() or 1)
        assert int(r["faces"]) == (int(r["divisions"])+1)**2
        for name, value in (("analytic", full), ("original", full), ("twoPoint", base),
                            ("limited", limited), ("naturalScale", scale)):
            equal(r[name], value, scale)
        equal(r["normalizedError"], abs(limited-full)/scale, D(1))
    failed, = result["failedCuts"]
    assert (failed["axis"], failed["field"], failed["shear"]) == (0, "y", .5)
    assert not failed["finestPass"] and not failed["orderPass"] and not result["passed"] and not receipt["methodGate"]
    # 在完整给定温差区间内另造 θ=0.15+0.1Y。端值一致而截面法向分量非零。
    reference = json.loads((ROOT/"cases/short_sphere/supportedFlux65.json").read_bytes())
    assert 0 < .15 < .25 < reference["frames"][-1]["maxScalar"]
    g = (D(0), D(".1"), D(0))
    cut_normal = (D(1), D("-.5"), D(0))
    _, limited_bounded, _ = projected((1,0,0), cut_normal, g, k)
    exact_bounded = -k*dot(cut_normal, g)
    assert exact_bounded == D(".05")*k and limited_bounded == 0
    for e in failed["errors"]:
        equal(e, 1/D(5).sqrt(), D(1))
    assert time.perf_counter()-started < 60
    return dict(cuts=len(cuts), starEdges=len(native), starRows=len(rows), worstEdge=worst_edge,
                residualErrorByField={key:float(v) for key,v in actual_ledger.items()},
                boundedAffineRange=[.15,.25], boundedAffineAnalyticCut=float(exact_bounded),
                boundedAffineLimitedCut=float(limited_bounded), finestAnalyticRelativeError=float(abs(limited_bounded/exact_bounded-1)),
                naturalScaleRelativeError=float(1/D(5).sqrt()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=EVIDENCE)
    parser.add_argument("--result", type=Path, default=ROOT/"cases/short_sphere/supportedDiffusion65.json")
    args = parser.parse_args()
    print(json.dumps(verify(args.result, args.raw), indent=2, allow_nan=False))
