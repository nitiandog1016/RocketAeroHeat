"""独立解包参考资料，核对工况记录、CSV 摘要、行序和六列原值。"""
import argparse
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
from zipfile import ZipFile

for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HEADER = b"sampleID,authorX_m,authorY_m,authorZ_m,authorPressure_Pa,authorHeatFlux_Wm2,authorWallTemperature_K"


def verify(bundle, direct):
    wanted = json.loads((ROOT / "cases/rfz_rec/deg190_source_record.json").read_bytes())
    with ZipFile(bundle) as z:
        assert sorted(z.namelist()) == ["samples.csv", "source.json"]
        meta = json.loads(z.read("source.json"))
        stored = meta.pop("exportedCSV")
        assert meta == wanted and stored["filename"] == "samples.csv"
        data = z.read(stored["filename"])
    digest = sha256(data).hexdigest()
    assert digest == stored["sha256"] and (direct is None or data == direct.read_bytes())
    with BytesIO(data) as stream:
        assert stream.readline().strip() == HEADER
        table = np.loadtxt(stream, delimiter=",", dtype="<f8")
    n = wanted["samples"]["count"]
    assert table.shape == (n, 7) and stored["rows"] == n and np.array_equal(table[:,0], np.arange(1,n+1))
    values = table[:,1:].tobytes(order="F")
    assert sha256(values).hexdigest() == wanted["samples"]["orderedValuesLEFloat64SHA256"]
    return dict(dataset=wanted["dataset"], samples=n, binary64ValuesVerified=n*6,
                directCSVMatched=direct is not None, csvSHA256=digest,
                readmeMach=wanted["readme"]["nominalMach"], workbookMach=wanted["workbook"]["mach"],
                paperCase=wanted["paper"]["caseNumber"], workbookCase=wanted["workbook"]["caseNumber"],
                pressureDatum=wanted["samples"]["pressureDatum"], heatFluxPositiveDirection=wanted["samples"]["heatFluxPositiveDirection"])


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("bundle",type=Path);p.add_argument("--csv",type=Path)
    args=p.parse_args()
    print(json.dumps(verify(args.bundle,args.csv),indent=2,allow_nan=False))
