"""核对独立安装的冷启动链、网格步骤和实际五步工作台所用文件及逐面结果。"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha
from thermistorBead import AREA as PAPER_AREA
from thermistorSeedRecovery import BASE_FRAME, COARSE_SEED, RELAX_FRAME
from thermistorStartupCache import BODY_SHA, COARSE_MESH, COARSE_NODES, GRID_SPEC, OBSERVER_SPEC


def read(file):
    return json.loads(file.read_bytes())


def save(file, value):
    with file.open("x",encoding="utf8") as f:
        json.dump(value,f,ensure_ascii=False,indent=2)


def preview(install):
    folder=install/"outputs";cache=folder/"thermistor_inputs/B6-56";view=read(folder/"preview.json")
    assert (cache/"initial").exists() is False and not (cache/"grid56/seed_00000.csv").exists()
    assert not any(p.name.startswith("v") for p in folder.iterdir())
    assert Path(view["directory"])==cache/"grid56" and view["bodySHA256"]==BODY_SHA
    assert view["meshSHA256"]==GRID_SPEC["gridSHA256"]==sha(cache/"grid56/mesh.su2")
    assert view["nodeMapSHA256"]==GRID_SPEC["nodeMapSHA256"]==sha(cache/"grid56/nodes.npz")
    v,f=np.array(view["preview"]["V"]),np.array(view["preview"]["F"])-1
    tri=v[f];assert len(f)==6048 and len(v)==3026 and np.all((tri.mean(1)*np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])).sum(1)>0)
    receipts={}
    for name in ("mesh48","grid56"):
        work=Path(read(cache/name/"installed.json")["sourceDirectory"])
        rec=read(work/"debugger_process.json");assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0
        assert not rec["stopReason"] and not rec["survivingProcessIDs"] and all(rec["inputsStable"].values())
        receipts[name]=dict(seconds=rec["elapsedSeconds"],peakRSSBytes=rec["peakFamilyRSSBytes"],receiptSHA256=sha(work/"debugger_process.json"))
    return dict(install=str(install),nativeSolves=0,previewVertices=len(v),previewTriangles=len(f),
                wallFaces=view["wallFaces"],cells=view["cells"],nodes=view["nodes"],meshSHA256=view["meshSHA256"],
                nodeMapSHA256=view["nodeMapSHA256"],seedPresent=False,receipts=receipts)


def cold(install):
    cache=install/"outputs/thermistor_inputs/B6-56"
    published=cache/"initial"
    meta=published/"installed.json"
    source=Path(read(meta)["sourceDirectory"]) if meta.exists() else published
    assert read(source/"seed_recovery.json")["source"]=="internal-cold"
    records={}
    for name in ("axis3","sector3","strong3","base6","relax6"):
        directory=source/name
        official=name in ("axis3","sector3","strong3")
        receipt=directory/("engine_process.json" if official else "debugger_process.json")
        rec=read(receipt)
        assert rec["status"]=="terminal" and rec["inputsStable"] and all(rec["inputsStable"].values())
        if official:
            assert rec["exitCode"]==0 and not rec["cancelled"] and not rec["timedOut"]
            assert rec["outputSHA256"]["restart.csv"]==sha(directory/"restart.csv")
        else:
            assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and not rec["survivingProcessIDs"]
            assert "Exit Success (SU2_CFD)" in (directory/"debugger.log").read_text(errors="replace")
        for file in rec["inputSHA256"]:
            Path(file).relative_to(install)
        records[name]=dict(frameSHA256=sha(directory/"restart.csv"),elapsedSeconds=rec["elapsedSeconds"],receiptSHA256=sha(receipt))
    assert records["base6"]["frameSHA256"]==BASE_FRAME and records["relax6"]["frameSHA256"]==RELAX_FRAME
    assert records["strong3"]["frameSHA256"]=="5f7eef6507a851ac0e6bc5c832dd2207a92bab55f0b83ea2ddfd4ace644957ba"
    assert sha(published/"seed_00000.csv")==COARSE_SEED
    assert sha(cache/"grid56/seed_00000.csv")==GRID_SPEC["seedSHA256"]
    assert sha(cache/"mesh48/nodes.npz")==COARSE_NODES and sha(cache/"mesh48/mesh.su2")==COARSE_MESH
    assert not any((install/"outputs"/n).exists() for n in ("v283_validation","v289_validation","v292_validation","v293_validation"))
    return dict(install=str(install),stages=records,coarseSeedSHA256=COARSE_SEED,fineSeedSHA256=GRID_SPEC["seedSHA256"],
                exactPhysicalStart=True)


def result(install, ui):
    u=read(ui);job=Path(u["directory"]);job.relative_to(install/"outputs/v294_validation/workbench")
    assert u["mode"]=="run" and u["webglError"]==0 and u["faces"]==5216 and Path(u["page"]).is_file()
    assert u["physicalTimeS"]==7.80936086544745e-6 and u["pickedFace"]==5211 and u["invalidatedOnGridChange"]
    # 将界面显示的终场绑定到同一次发起计算的监督回执及该任务的实际输入。
    outer_path=ui.parent/"debugger_process.json"
    outer=read(outer_path)
    assert outer["status"]=="terminal" and outer["debuggerExitCode"]==0
    assert not outer["stopReason"] and not outer["supervisorError"] and not outer["survivingProcessIDs"]
    assert outer["inputsStable"] and all(outer["inputsStable"].values())
    request=read(job/"request.json")
    assert Path(request["directory"])==job and request["seedSHA256"]==GRID_SPEC["seedSHA256"]
    assert request["meshSHA256"]==GRID_SPEC["gridSHA256"]
    cache=install/"outputs/thermistor_inputs/B6-56"
    initial=cache/"initial/seed_00000.csv"
    assert read(initial.parent/"installed.json")["source"]=="verified-internal-checkpoint"
    assert sha(initial)==COARSE_SEED
    initial_fine=cache/"grid56/seed_00000.csv"
    assert sha(initial_fine)==sha(job/"first01/seed_00000.csv")==request["seedSHA256"]
    prolong_receipts=[p for p in cache.glob("prolong-*/debugger_process.json")
                      if (p.parent/initial.name).is_file() and initial_fine.samefile(p.parent/initial.name)]
    assert len(prolong_receipts)==1
    prolong=read(prolong_receipts[0])
    assert prolong["status"]=="terminal" and prolong["debuggerExitCode"]==0 and not prolong["stopReason"]
    assert not prolong["supervisorError"] and not prolong["survivingProcessIDs"] and all(prolong["inputsStable"].values())
    assert prolong["inputSHA256"][str(initial)]==COARSE_SEED
    costs={}
    for p in sorted(job.glob("*/debugger_process.json")):
        rec=read(p);assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and not rec["stopReason"]
        assert not rec["supervisorError"] and not rec["survivingProcessIDs"] and all(rec["inputsStable"].values())
        for name,digest in rec["inputSHA256"].items():
            path=Path(name);path.relative_to(install)
            assert path.is_file() and sha(path)==digest,path
        costs[p.parent.name]=dict(elapsedSeconds=rec["elapsedSeconds"],peakFamilyRSSBytes=rec["peakFamilyRSSBytes"],receiptSHA256=sha(p))
    progress=read(job/"progress.json")
    assert len(costs)==11 and progress["status"]=="accepted"
    assert read(Path(u["resultReceipt"]))==progress["result"]
    bundle=install/"cases/thermistor_b/startup_6ms_nearwall56.zip"
    with zipfile.ZipFile(bundle) as z:
        raw=z.read("wall.csv")
        reference=json.loads(z.read("summary.json"))
    wall=job/"second_energy/wall.csv"
    assert raw==wall.read_bytes() and hashlib.sha256(raw).hexdigest()==u["wallSHA256"]
    assert sha(job/"second01/restart_00002.csv")==u["sourceFrame"]
    assert not any((install/"outputs"/n).exists() for n in ("v283_validation","v289_validation","v292_validation","v293_validation"))
    observer=install/"outputs/thermistor_inputs/B6-56/observer"
    b=read(observer/"build.json");assert b["origin"]=="rebuilt" and b["verifiedInputs"]==302
    assert sha(observer/"observe_bead.exe")==OBSERVER_SPEC["observerSHA256"]
    src=Path(b["sourceDirectory"]);compile_rec,link_rec=[read(src/(s+"_process.json")) for s in ("compile","link")]
    assert all(r["exitCode"]==0 and not r["stopReason"] and not r["livePids"] for r in (compile_rec,link_rec))
    comparison=read(job/"compare/summary.json")
    surface=read(job/"second_energy/accepted.json")["surface"]
    energy=surface["energy"]
    assert comparison["withinComparisonThresholds"] and comparison["fineSurface"]==surface
    assert surface["fullSurfaceComputed"] and surface["areaM2"]==reference["areaM2"]
    assert abs(u["areaM2"]-surface["areaM2"])<1e-15 and abs(u["intoWallW"]+energy["wallReleaseW"])<1e-13
    assert surface["peakOutwardFaceID"]+1==u["pickedFace"]
    first=read(job/"first_state/verified.json")
    second=read(job/"second_state/verified.json")
    assert first["passed"] and second["passed"] and first["frameSHA256"]==second["oldFrameSHA256"]
    assert first["oldFrameSHA256"]==request["seedSHA256"] and second["frameSHA256"]==u["sourceFrame"]
    # 用论文名义面积换算整珠系数，同时保留用于原面计算覆盖率的三角面积。
    condition=request["condition"]
    measured=read(HERE.parent/"v287/evidence/measurements.json")["heatExperiment"]
    assert measured["speedMps"]==condition["speedMps"] and measured["measuredQuantity"]=="whole-bead mean h"
    mean_h=energy["wallReleaseW"]/(PAPER_AREA*(condition["wallK"]-condition["airK"]))
    experiment=dict(doi=measured["doi"],figure=measured["figure"],printedPage=measured["printedPage"],
                    nominalAreaM2=PAPER_AREA,startupMeanHWm2K=mean_h,measuredMeanHWm2K=measured["hWm2K"],
                    measuredDeltaK=measured["deltaK"],graphReadRangeWm2K=measured["graphReadRangeWm2K"],
                    reportedPointUncertainty=measured["experimentalPointUncertainty"],
                    relativeDifference=mean_h/measured["hWm2K"]-1)
    return dict(workbench=u,condition=request["condition"],acceptedStages=costs,wallByteIdenticalToReference=True,
                wholeBeadComparison=experiment,
                pressurePa=surface["wallPressurePa"],heatFluxIntoWallWm2=surface["heatFluxWm2"],
                nearWallFluidK=surface["nearWallK"],jointExperimentVerifiedAreaM2=surface["jointExperimentVerifiedAreaM2"],
                trailingPassingRows=[first["trailingPassingRows"],second["trailingPassingRows"]],
                coarseSeedSHA256=COARSE_SEED,fineSeedSHA256=request["seedSHA256"],prolongSeconds=prolong["elapsedSeconds"],
                matlabElapsedIncludingSolvesS=outer["elapsedSeconds"],matlabFamilyPeakRSSBytes=outer["peakFamilyRSSBytes"],
                matlabReceiptSHA256=sha(outer_path),
                observerSHA256=b["observerSHA256"],observerSourceSHA256=b["sourceSHA256"],
                compileSeconds=compile_rec["elapsedSeconds"],linkSeconds=link_rec["elapsedSeconds"],
                wallReleaseW=energy["wallReleaseW"],imbalanceW=energy["imbalanceW"])


def archive_evidence(target, install, ui):
    """只归档本次小型原始回执、诊断与实测源码差异，大型体网和场由指纹绑定。"""
    job = Path(read(ui)["directory"])
    area = install.parent
    cache = install / "outputs/thermistor_inputs/B6-56"
    cold_install = Path(read(HERE / "evidence/cold.json")["install"])
    cold_cache = cold_install / "outputs/thermistor_inputs/B6-56"
    initial = cold_cache / "initial"
    if (initial / "installed.json").is_file():
        initial = Path(read(initial / "installed.json")["sourceDirectory"])
    sources = [ui, ui.parent / "debugger_process.json", ui.parent / "debugger.log",
               job / "request.json", job / "progress.json",
               area / "cold02_executed_cache.py", area / "cold.log",
               initial / "seed_recovery.json"]
    for name in ("first", "second"):
        sources.extend([job / (name + "_state/verified.json"), job / (name + "_energy/accepted.json")])
    sources.extend(job.glob("*/debugger_process.json"))
    for name in ("axis3", "sector3", "strong3", "base6", "relax6"):
        sources.extend((initial / name).glob("*_process.json"))
    for name in ("mesh48", "grid56", "observer"):
        installed = cache / name / ("build.json" if name == "observer" else "installed.json")
        sources.append(installed)
        sources.extend(Path(read(installed)["sourceDirectory"]).glob("*_process.json"))
    sources.extend(cache.glob("prolong-*/debugger_process.json"))
    sources.extend(cold_cache.glob("initial-*/axis3/supervisor.log"))
    for folder in ("mesh_ui01", "mesh_ui02"):
        sources.extend((area / folder).glob("debugger*"))
        sources.extend((area / folder).glob("*.m"))
    for name in ("scripts/thermistorStartupCache.py", "scripts/thermistorStartupObserver.py",
                 "scripts/runThermistorStartupJob.py", "prepareThermistorWorkbenchMesh.m",
                 "rocketAeroHeatWorkbench.m"):
        sources.append(install / name)
    with zipfile.ZipFile(target, "x", compression=zipfile.ZIP_DEFLATED) as bundle:
        for file in sorted(set(sources)):
            relative = file.resolve().relative_to(area)
            assert file.stat().st_size < 1024**2, file
            bundle.write(file, relative.as_posix())
        # 运行后更新安装清单时，仍可从本次原始回执和快照核对测试驱动。
        driver = ui.parent / "executed_install_check.py"
        expected_driver = read(ui.parent / "debugger_process.json")["inputSHA256"][str(HERE / "install_check.py")]
        assert sha(driver) == expected_driver
        bundle.write(driver, driver.relative_to(area).as_posix())


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("mode",choices=("preview","cold","result"))
    cli.add_argument("install",type=Path)
    cli.add_argument("out",type=Path)
    cli.add_argument("--ui",type=Path)
    cli.add_argument("--archive",type=Path,help="同时保存此次实际运行的小型原始回执")
    a=cli.parse_args();inst=a.install.resolve();inst.relative_to(ROOT/"outputs/v296_validation")
    if a.mode=="result" and a.ui is None:cli.error("端到端核对需要 --ui 指定五步工作台回执")
    if a.archive and a.mode!="result":cli.error("完整回执归档用于端到端验收")
    if a.mode=="preview":r=preview(inst)
    elif a.mode=="cold":r=cold(inst)
    else:
        r=result(inst,a.ui)
        if a.archive:archive_evidence(a.archive,inst,a.ui.resolve())
    save(a.out,r)
