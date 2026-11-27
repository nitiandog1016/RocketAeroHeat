"""对照原面归档，核对工作台实际导出的全表面、三条截面曲线与粗细网来源。"""
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

import numpy as np
from refined import DOC, ROOT, read, save, sha
sys.path.insert(0, str(DOC.parent/"v351"))
from verify_scopes import EXPORT


def check(ui, web):
    result = read(ui/"ui.json")
    case = ROOT/"cases/simple_launcher/wake"
    seal = read(case/"receipt.json")
    assert sha(case/"wind_tunnel.zip") == seal["sha256"]
    assert sha(case/"wake_heat_grid.zip") == seal["wakeHeatGridSHA256"]
    with ZipFile(case/"wind_tunnel.zip") as z:
        coarse = np.loadtxt(z.open("wall.csv"),delimiter=",",skiprows=1)
    with ZipFile(case/"wake_heat_grid.zip") as z:
        fine = np.loadtxt(z.open("fine_wall.csv"),delimiter=",",skiprows=1)
        pair = json.loads(z.read("comparison.json"))
    time = pair["condition"]["elapsedPhysicalTimeS"]
    assert all(result[k] == pair[k] for k in ("coarseWallSHA256","fineWallSHA256","flowSHA256"))
    np.testing.assert_array_equal(coarse[:,:3],fine[:,:3])
    values = dict(heatFlux=coarse[:,3],wakeFineHeat=fine[:,3],wakeHeatDifference=fine[:,3]-coarse[:,3])
    # 可指定发布版网页脚本目录；独立重做几何截面，核对工作台导出的同一批原面。
    cut_program = """
const fs=require('node:fs'),path=require('node:path');
const [folder,lib]=process.argv.slice(1);
const raw=fs.readFileSync(path.join(folder,'case-data.js'),'utf8').trim();
const data=JSON.parse(raw.slice('window.ROCKET_CASE='.length,-1));
const tool=require(path.join(lib,'section.js'));
const cut=tool.geometry(data.vertices,data.faces,0.7);
const result={};
for(const name of ['heatFlux','wakeFineHeat','wakeHeatDifference']){
  const field=data.fields.find(f=>f.key===name);
  result[name]=tool.sample(cut,field).map(s=>[s.face+1,s.value]);
}
process.stdout.write(JSON.stringify(result));
"""
    section = subprocess.run(["node","-e",cut_program,str((ui/"workbench").resolve()),str(web.resolve())],
                             capture_output=True,check=True,timeout=45)
    samples = json.loads(section.stdout)
    counts = {}
    for key, field in values.items():
        proc = subprocess.run(["node","-e",EXPORT,str((ui/"workbench").resolve()),str(web.resolve()),key],
                              stdout=subprocess.PIPE,check=True,timeout=45)
        rows = list(csv.DictReader(io.StringIO(proc.stdout.decode("utf-8-sig"))))
        assert len(rows) == 9216
        assert all(row["单位"] == "W/m²" and row["算法"] and row["场量标识"] == key and
                   float(row["时间 / s"]) == time for row in rows)
        np.testing.assert_array_equal([int(next(iter(r.values()))) for r in rows],fine[:,0])
        np.testing.assert_allclose([float(r["数值"]) for r in rows],field,rtol=0,atol=1e-13)
        with (ui/(key+"_section.csv")).open(encoding="utf-8-sig",newline="") as f:
            cut = list(csv.DictReader(f))
        assert cut and all(float(r["elapsed_physical_time_s"]) == time and r["field_key"] == key for r in cut)
        np.testing.assert_allclose([float(r["value"]) for r in cut],field[[int(r["face_id"])-1 for r in cut]],rtol=0,atol=1e-13)
        actual = sorted((int(r["face_id"]),float(r["value"])) for r in cut)
        np.testing.assert_allclose(sorted(samples[key]),actual,rtol=0,atol=1e-13)
        counts[key] = dict(fullRows=len(rows),sectionRows=len(cut))
    with (ui/"wall_pair.csv").open(encoding="utf-8-sig",newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows)==len(fine) and all(r["joint_verified"]=="0" and r["source_traceable"]=="1" for r in rows)
    np.testing.assert_array_equal([int(r["face_id"]) for r in rows],coarse[:,0])
    np.testing.assert_allclose([float(r["pressure"])*1000 for r in rows],coarse[:,2],atol=1e-8,rtol=0)
    np.testing.assert_allclose([float(r["heat_flux"]) for r in rows],coarse[:,3],atol=1e-13,rtol=0)
    return dict(ui=result,exports=counts,pairedRows=len(rows),jointExperimentalAreaM2=0,
                exportJSSHA256=sha(web/"export.js"),sectionJSSHA256=sha(web/"section.js"))


if __name__ == "__main__":
    ui, web, output = map(Path,sys.argv[1:])
    save(output,check(ui,web))
