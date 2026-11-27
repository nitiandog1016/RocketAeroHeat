"""独立回读普通工作台的全周复核 CSV、截面及联合面表，核对发布包。"""
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(DOC.parent/"v351"))
from verify_scopes import EXPORT


def check(folder,web=None):
    ui=json.loads((folder/"ui.json").read_bytes())
    job=Path(ui["runDirectory"])
    source=ROOT/"cases/simple_launcher" if ui["mode"]=="archive" else job/"result"
    pack=source/"angular_heat.zip"
    assert hashlib.sha256(pack.read_bytes()).hexdigest()==json.loads((source/"receipt.json").read_bytes())["angularHeatSHA256"]
    with ZipFile(pack) as z:
        meta=json.loads(z.read("comparison.json"))
        raw=z.read("angular_wall.csv")
    assert hashlib.sha256(raw).hexdigest()==meta["csvSHA256"]==ui["angularCSV"]
    assert meta["fieldSHA256"]==ui["angularFieldSHA256"] and meta["flowSHA256"]==ui["flowSHA256"]
    original=np.loadtxt(io.BytesIO(raw),delimiter=",",skiprows=1)
    np.testing.assert_array_equal(original[:,0],np.arange(1,9217))
    page=(folder/"workbench").resolve()
    run=subprocess.run(["node","-e",EXPORT,str(page),str(web or page),"launcherAngularHeat"],capture_output=True,check=True,timeout=45)
    reader=csv.DictReader(io.StringIO(run.stdout.decode("utf-8-sig")))
    values=[]
    for face,row in enumerate(reader,1):
        assert int(row[reader.fieldnames[0]])==face and row["单位"]=="W/m²"
        assert "48 方位" in row["算法"] and meta["fieldSHA256"] in row["逐面元适用范围"]
        values.append(float(row["数值"]))
    np.testing.assert_allclose(values,original[:,2],rtol=0,atol=1e-14)
    with (folder/"section.csv").open(encoding="utf-8-sig",newline="") as stream:
        cut=list(csv.DictReader(stream))
    np.testing.assert_allclose([float(r["value"]) for r in cut],original[[int(r["face_id"])-1 for r in cut],2],rtol=0,atol=1e-13)
    assert len(cut)==ui["sectionRows"] and abs(original[:,1]@values+meta["outwardW"])<1e-12
    spec=importlib.util.spec_from_file_location("v350_verify",DOC.parent/"v350/verify_delivery.py")
    validator=importlib.util.module_from_spec(spec);spec.loader.exec_module(validator)
    primary=validator.verify(folder,True)
    return dict(ui=ui,fullSurfaceCSVRows=len(values),sectionRows=len(cut),outwardW=meta["outwardW"],primary=primary)


if __name__=="__main__":
    print(json.dumps(check(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve() if len(sys.argv)>2 else None),indent=2))
