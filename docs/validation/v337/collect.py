"""从已验收的三个稳态和两份同刻保存帧提取七站截面，供工作台归档。"""
import argparse
import json
import os
from pathlib import Path
import sys

for var in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[var]="1"
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from shortSupportStationGradient import measurement_grid,reconstruct
from shortSphere65 import read,sha

CASE=ROOT/"cases/short_sphere"
VOLUME=ROOT/"outputs"


def collect(out):
    normal=read(CASE/"supportedNearWall65.json")
    angular=read(CASE/"supportedAngular65.json")
    pair=read(CASE/"supportedTime65.json")
    rows=[]
    app=VOLUME/"v336_validation/clean01"
    job=app/"outputs/short65_time_6dea1a68-073d-4129-867d-3772b804b131"
    origin=read(job/"construction.json")
    assert sha(job/"construction.json")==read(job/"completed.json")["constructionSHA256"]
    for label,prep,steady,archive,expected,frame_key in (
        ("steady5",VOLUME/"v325_validation/prep02/prepared",VOLUME/"v329_validation/run01","supportedSurface65.zip",normal["coarseArchiveSHA256"],normal["coarseFrame"]),
        ("steady10",VOLUME/"v330_validation/prepared01",VOLUME/"v330_validation/run01","supportedFine65.zip",normal["fineArchiveSHA256"],normal["fineFrame"]),
        ("steady128",app/origin["prepared"],app/origin["steady"],"supportedAngularSurface65.zip",angular["fineArchiveSHA256"],angular["fineFrame"])):
        assert sha(CASE/archive)==expected
        grid=measurement_grid(prep)
        end=read(steady/"converged.json")
        assert read(steady/"thermal.json")["passed"] and end["frame"]==frame_key
        stage=steady/end["stage"]
        terminal=read(stage/"debugger_process.json")
        assert terminal["status"]=="terminal" and terminal["debuggerExitCode"]==0 and all(terminal["inputsStable"].values())
        result=reconstruct(grid,stage/"restart_000199.csv",CASE/archive)
        rows.append(dict(label=label,**result))
    for branch,step in (("full",1),("half",2)):
        stage=job/"run"/branch/f"step{step}"
        assert read(stage/"native.json")["frame"]==pair[branch+"Frame"]
        assert abs(read(stage/"energy.json")["relativeBalance"])<.005
        data=job/"export"/branch/"surface.zip"
        assert sha(data)==read(job/"export/pair.json")[branch+"Archive"]
        result=reconstruct(grid,stage/f"restart_{step:05d}.csv",data)
        rows.append(dict(label=branch,**result))
    with out.open("x",encoding="utf8",newline="\n") as f:
        json.dump(dict(schema="short65-station-gradient-v1",sets=rows),f,ensure_ascii=False,allow_nan=False,indent=1)
    print([(g["label"],sum(row["status"]=="available" for row in g["rows"])) for g in rows])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("output",type=Path);collect(p.parse_args().output)
