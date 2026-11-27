"""将完整面表、工作台的逐变量导出和截面逐面比较，检查每个数值的时刻与来源。"""
import csv
import hashlib
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


def check(ui,source,web):
    result=json.loads((ui/"ui.json").read_bytes())
    pack=source/"wind_tunnel.zip"
    seal=json.loads((source/"receipt.json").read_bytes())
    assert hashlib.sha256(pack.read_bytes()).hexdigest()==seal["sha256"]
    with ZipFile(pack) as z:
        summary=json.loads(z.read("summary.json"));raw=z.read("wall.csv")
    condition=summary["condition"];time=condition["elapsedPhysicalTimeS"]
    assert hashlib.sha256(raw).hexdigest()==result["wallSHA256"]==summary["wallSHA256"]
    assert result["sourceFlow"]==summary["sourceFrameSHA256"] and result["sourceHeat"]==summary["thermalFrameSHA256"]
    native=np.loadtxt(io.BytesIO(raw),delimiter=",",skiprows=1)
    count=len(native);np.testing.assert_array_equal(native[:,0],np.arange(1,9217))
    with (ui/"wall_pair.csv").open(encoding="utf-8-sig",newline="") as f:paired=list(csv.DictReader(f))
    assert len(paired)==count and all(int(r["source_traceable"])==1 and int(r["joint_verified"])==0 for r in paired)
    np.testing.assert_array_equal([int(r["face_id"]) for r in paired],native[:,0])
    np.testing.assert_allclose([float(r["pressure"])*1000 for r in paired],native[:,2],rtol=0,atol=1e-8)
    np.testing.assert_allclose([float(r["heat_flux"]) for r in paired],native[:,3],rtol=0,atol=1e-13)
    page=(ui/"workbench").resolve()
    quantities={"pressure":(native[:,2]/1000,"kPa"),"pressureGauge":(native[:,2]-101325,"Pa"),
            "heatFlux":(native[:,3],"W/m²"),"nearWallCellTemperature":(native[:,4],"K"),
            "wallInput":(np.full(count,300.001),"K")}
    for key,(values,unit) in quantities.items():
        run=subprocess.run(["node","-e",EXPORT,str(page),str(web.resolve()),key],check=True,capture_output=True,timeout=45)
        table=list(csv.DictReader(io.StringIO(run.stdout.decode("utf-8-sig"))))
        assert len(table)==count and all(r["单位"]==unit and abs(float(r["时间 / s"])-time)<1e-12 for r in table)
        assert all("0.600 s" in r["工况"] and "稳态轴对称" not in r["工况"] for r in table)
        np.testing.assert_allclose([float(r["数值"]) for r in table],values,atol=1e-12,rtol=0)
    sections={}
    for key in ("heatFlux","pressureGauge"):
        with (ui/(key+"_section.csv")).open(encoding="utf-8-sig",newline="") as f:cut=list(csv.DictReader(f))
        assert cut and all(abs(float(r["elapsed_physical_time_s"])-time)<1e-12 for r in cut)
        idx=[int(r["face_id"])-1 for r in cut]
        np.testing.assert_allclose([float(r["value"]) for r in cut],quantities[key][0][idx],rtol=0,atol=1e-12)
        sections[key]=len(cut)
    return dict(result=result,rowsPerField=count,fields=list(quantities),sectionRows=sections,
                pairedRows=len(paired),coveredAreaM2=float(native[:,1].sum()),outwardW=float(-native[:,1]@native[:,3]))


if __name__=="__main__":
    print(json.dumps(check(*(Path(x).resolve() for x in sys.argv[1:])),indent=2))
