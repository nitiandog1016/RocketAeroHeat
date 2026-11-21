"""把工作台导出的实际离线页面与同次热敏珠计算的原面表独立逐行比对。"""
import csv
import hashlib
import json
from pathlib import Path
import sys


def audit(receipt):
    ui = json.loads(receipt.read_bytes())
    page = Path(ui["page"])
    script = page.parent / "case-data.js"
    raw = script.read_text("utf8").strip()
    assert page.is_file() and raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    data = json.loads(raw[len("window.ROCKET_CASE="):-1])
    identity = data["evidenceIdentity"]
    assert identity["runDirectory"] == ui["directory"] and identity["startupAction"] == ui["mode"]
    assert identity["sourceFrame"]["sha256"] == ui["sourceFrame"]
    wall = Path(ui["directory"]) / "second_energy/wall.csv"
    assert hashlib.sha256(wall.read_bytes()).hexdigest() == identity["sourceFrame"]["wallSHA256"] == ui["wallSHA256"]
    with wall.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))
    n = len(rows)
    assert n == 5216 and [int(r["faceID"]) for r in rows] == list(range(n))
    fields = {f["key"]: f for f in data["fields"]}
    errors = {}
    for key, column, scale in (("pressure","pressurePa",1000), ("pressureGauge","pressureGaugePa",1),
                               ("heatFlux","heatFluxWm2",1000), ("nearWallCellTemperature","nearWallCellTemperatureK",1)):
        field = fields[key]
        assert len(field["values"]) == n and field["status"] == [1]*n
        errors[key] = max(abs(display*scale - float(r[column])) for display, r in zip(field["values"],rows))
        assert errors[key] < 1e-9, (key, errors[key])
    assert fields["wallInput"]["values"] == [303.15]*n and fields["wallInput"]["status"] == [2]*n
    for key in ("wallSolved", "edgeTemperature", "recoveryTemperature"):
        assert fields[key]["values"] == [None]*n and fields[key]["status"] == [0]*n
    assert data["resultCondition"]["elapsedPhysicalTimeS"] == ui["physicalTimeS"]
    return dict(page=str(page), caseDataSHA256=hashlib.sha256(script.read_bytes()).hexdigest(),
                wallSHA256=ui["wallSHA256"], checkedFaces=n, maxErrorInSI=errors, givenWallK=303.15,
                unknownFields=["wallSolved","edgeTemperature","recoveryTemperature"])


if __name__ == "__main__":
    result = audit(Path(sys.argv[1]))
    with Path(sys.argv[2]).open("x", encoding="utf8") as out:
        json.dump(result, out, ensure_ascii=False, indent=2)
