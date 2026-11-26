"""从作者发布的 NPY 独立逐行复核参考 CSV，包括负值、坐标、样点顺序及双精度。"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "cases/rfz_rec"
FIELDS = ("X", "Y", "Z", "p", "q", "T")
HEADER = "sampleID,authorX_m,authorY_m,authorZ_m,authorPressure_Pa,authorHeatFlux_Wm2,authorWallTemperature_K"


def verify(csv_file, source):
    coords = json.loads((ROOT / "docs/validation/v342/evidence/reference.json").read_bytes())
    boundary = json.loads((CASE / "deg190_boundary.json").read_bytes())
    case = json.loads((CASE / "deg190_reference.json").read_bytes())
    hashes = dict(coords["coordinateSHA256"], **boundary["fieldSHA256"])
    if source is None:
        base = ROOT / "outputs/rfz_reference"
        source = base / "author_data_20261124/files/Aerodynamics/Stage_1" / case["dataset"]
        if not source.is_dir():
            source = base / "published_samples" / case["sourceCommit"]
    columns = []
    for key in FIELDS:
        file = source / (key+"_data.npy")
        assert sha256(file.read_bytes()).hexdigest() == hashes[key]
        values = np.load(file, allow_pickle=False)
        assert values.shape == (72706,) and values.dtype == "<f8" and np.isfinite(values).all()
        columns.append(values)
    reference = np.column_stack(columns)
    with csv_file.open(encoding="utf-8-sig") as stream:
        assert stream.readline().strip() == HEADER
        data = np.loadtxt(stream, delimiter=",")
    assert data.shape == (len(reference), 7)
    assert np.array_equal(data[:, 0], np.arange(1, len(data)+1))
    assert np.array_equal(data[:, 1:].view(np.uint64), reference.view(np.uint64))
    xyz, p, q, t = reference[:, :3], reference[:, 3], reference[:, 4], reference[:, 5]
    return dict(sourceCommit=case["sourceCommit"], dataset=case["dataset"], sourceSHA256=hashes,
                csvSHA256=sha256(csv_file.read_bytes()).hexdigest(), samples=len(data),
                distinctSourcePositions=len(np.unique(xyz, axis=0)), binary64ValuesMatched=int(reference.size),
                pressureRangePa=[p.min(), p.max()], heatFluxRangeWm2=[q.min(), q.max()],
                negativeHeatFluxSamples=int(np.sum(q < 0)), zeroHeatFluxSamples=int(np.sum(q == 0)),
                temperatureRangeK=[t.min(), t.max()], projectComputedFaces=0, projectComputedAreaM2=0.)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--source", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.csv_file, args.source), indent=2, allow_nan=False))
