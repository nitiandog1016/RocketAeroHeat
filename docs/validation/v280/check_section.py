"""用发布版网页脚本导出中跨截面，核对来源、原面编号与未收敛量的单位。"""
import csv
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from compare_step import HERE, ROOT, save


def check(page, out):
    with (HERE / "evidence/trial_faces.csv").open(encoding="utf-8-sig") as stream:
        faces = list(csv.DictReader(stream))
    frame = json.loads((HERE / "evidence/trial.json").read_bytes())["sourceFrameSHA256"]
    out.mkdir(parents=True, exist_ok=False)
    result = {}
    code = """
const path = require('path');
const [page, mods, key] = process.argv.slice(1);
global.window = {};
require(path.join(page,'case-data.js'));
const model = window.ROCKET_CASE, field = model.fields.find(f => f.key === key);
const section = require(path.join(mods,'section.js')), viewer = require(path.join(mods,'viewer.js'));
const cut = section.geometry(model.vertices,model.faces,0);
if(cut.components.length!==1 || !cut.components[0].closed) throw Error('open section');
process.stdout.write(section.csv(cut,field,'冻结工况',viewer.caseCaption(model,field,0),{model}));
"""
    with tempfile.TemporaryDirectory() as temp:
        release = Path(temp)
        for name in ("section.js", "viewer.js"):
            (release/name).write_bytes(subprocess.check_output(["git", "show", "v279:web/"+name], cwd=ROOT))
        for key, col, scale in (("trialPressure","trialPressurePa",1000), ("trialHeatFlux","trialHeatFluxWm2",1000)):
            for mods in (page, release):
                raw = subprocess.check_output(["node", "-e", code, str(page.resolve()), str(mods.resolve()), key])
                rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
                assert len(rows) == len({r["面元编号"] for r in rows}) == 720
                for row in rows:
                    face = faces[int(row["面元编号"])-1]
                    assert row["面元编号"] == face["faceID"] and row["状态"] == "未收敛试算"
                    assert frame in row["适用范围"] and "720 次" in row["适用范围"] and "自适应 CFL" in row["适用范围"]
                    assert row["算法"] and "303.15" in row["工况"]
                    field_col = next(k for k in row if k.startswith("未收敛试算"))
                    assert abs(float(row[field_col])*scale - float(face[col])) < 1e-8
            with gzip.GzipFile(filename=str(out / (key+".csv.gz")), mode="xb", mtime=0) as archive:
                archive.write(raw)
            result[key] = len(rows)
    save(out / "sections.json", dict(checkedRows=result, sourceFrame=frame, modules=["locally exported", "v279 release"]))


if __name__ == "__main__":
    check(Path(sys.argv[1]), Path(sys.argv[2]))
