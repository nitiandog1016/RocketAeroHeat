"""在当前安装内从 14D 银球的几何与温界建立冷态、验收基流并推进实际物理时间。"""
from pathlib import Path
import sys
import tempfile

from shortSphere65 import CASE, ROOT, read, sha, write_json
from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import exclusive
from installShort65Kernel import HOME_DIR, install, invoke
from buildShort65Programs import build, check_programs
from build_time import PREP_KEY
from run_shared import solve
from map_surface import mapped
from physical_step import setup, run_steps, inspect_step
from time_surface import export
from energy import endpoint

DATA=ROOT/"cases/short_sphere"
CACHE=ROOT/"outputs/short65_inputs/14D"
PLAN=ROOT/"docs/validation/v336/PLAN.md"
MESH="a56a4c3dee8fd8bb74bc23ae4e85b71eaf7d5552dd99fccf8a6c8cb6a40789cd"
COLD="82eb4bee259a454cbfaf3a3728b396cf30a92456b636a3a245c1f0f5f0e084c0"
TRANSFER="63c49c77d06cd2c1814d97f68769c62c3401d52748ad9863f0382fb8408a8609"


def directory(parent,name):
    path=(parent/name).resolve();path.relative_to(parent.resolve());return path


def phase(job,code,message):
    atomic_write_json(job/"progress.json",dict(stage=code,message=message));print(message,flush=True)


def prepare_grid(prep):
    """在独立进程内生成体网与冷态，限制首次建网的时间及进程族内存。"""
    import installGmshRuntime
    installGmshRuntime.main()
    invoke("mesh", [sys.executable, ROOT/"scripts/shortSupportMesh.py", prep, "--split-near", "--split-sectors"],
           ROOT, prep.parent, 180, 2, bin_directory=Path(sys.executable).parent)


def inputs(job):
    CACHE.mkdir(parents=True,exist_ok=True)
    with exclusive(CACHE):
        index=CACHE/"inputs.json";restore=index.is_file()
        if restore:
            keys=read(index);prep=directory(CACHE,keys["directory"])
            assert sha(prep/"prepared.json")==keys["prepared"]
        else:
            phase(job,"mesh","首次准备：生成十层近壁体网、五区交叠、给定温界和均匀冷初场")
            prep=Path(tempfile.mkdtemp(prefix="grid_",dir=str(CACHE)))/"prepared"
            prepare_grid(prep)
        spec=read(prep/"prepared.json")
        assert sha(prep/"prepared.json")==PREP_KEY
        assert spec["files"]["mesh.su2"]==MESH and spec["files"]["initial.csv"]==COLD and spec["caseSHA256"]==sha(CASE)
        assert all(sha(prep/n)==h for n,h in spec["files"].items())
        assert spec["supportedSHA256"]==sha(DATA/"supported65.zip") and sha(DATA/"supportedAngularTransfer65.npz")==TRANSFER
        if not restore:atomic_write_json(index,dict(directory=prep.relative_to(CACHE).as_posix(),prepared=sha(prep/"prepared.json")))
        programs=CACHE/"programs.json"
        if programs.is_file():
            entry=read(programs);native=directory(CACHE,entry["directory"])
            for k in ("steady","time"):assert sha(native/k/"build.json")==entry[k]
            check_programs(native,prep)
        else:
            phase(job,"programs","首次准备：构建并核对银球稳态、BDF1 及独立壁反力程序")
            if not (HOME_DIR/"kernel.json").is_file():install()
            native=Path(tempfile.mkdtemp(prefix="native_",dir=str(CACHE)))/"programs"
            build(native,prep)
            check_programs(native,prep)
            atomic_write_json(programs,dict(directory=native.relative_to(CACHE).as_posix(),steady=sha(native/"steady/build.json"),time=sha(native/"time/build.json")))
        base=CACHE/"base.json"
        if base.is_file():
            entry=read(base);steady=directory(CACHE,entry["directory"])
            assert sha(steady/"thermal.json")==entry["thermal"] and sha(steady.parent/"surface/surface.zip")==entry["surface"]
        else:
            phase(job,"steady","首次准备：从本安装冷态求解稳态基流，核对末窗、热稳定及壁面收支")
            steady=Path(tempfile.mkdtemp(prefix="base_",dir=str(CACHE)))/"steady"
            solve(steady,native/"steady",prepared=prep,plan=PLAN,release_driver=native/"steady/wire.exe")
            mapped(steady,steady.parent/"surface",prepared=prep,expected_mesh=MESH,transfer_file=DATA/"supportedAngularTransfer65.npz",transfer_sha=TRANSFER)
        stages=check_steady(steady,prep,native)
        if not base.is_file():
            atomic_write_json(base,dict(directory=steady.relative_to(CACHE).as_posix(),thermal=sha(steady/"thermal.json"),
                                      surface=sha(steady.parent/"surface/surface.zip")))
        return prep,native,steady,stages


def check_steady(steady, prepared, native):
    """复用稳态前将实际端点的输入、程序、收支和保存帧连接到当前安装。"""
    stages = read(steady / "converged.json")
    last = directory(steady, stages["stage"])
    record = read(last / "record.json")
    proc = read(last / "debugger_process.json")
    exe = native / "steady/shared.exe"
    thermal = read(steady / "thermal.json")
    assert thermal["passed"] and thermal["early"]["sourceFrame"] == stages["frame"]
    assert thermal["caseSHA256"] == sha(CASE) and thermal["meshSHA256"] == MESH
    assert proc["status"] == "terminal" and proc["debuggerExitCode"] == 0
    assert not (proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
    assert all(proc["inputsStable"].values()) and proc["command"][0] == str(exe)
    for file in (exe, CASE, prepared / "prepared.json", last / "case.cfg", last / "wire_boundary.dat",
                 last / "mesh.su2", last / "initial_flow.csv", Path(solve.__code__.co_filename)):
        assert sha(file) == proc["inputSHA256"][str(file)]
    assert record["exe"] == sha(exe) and max(record["window"]) < -8 and record["rawMax20"] < 1e-8
    assert sha(last / "restart_000199.csv") == stages["frame"] == record["end"]
    return stages


def run(job,completed=None):
    if completed is None:
        prep,native,steady,stages=inputs(job);base=stages["frame"];work=job/"run"
        origin=dict(schema="short65-cold-start-v1",caseSHA256=sha(CASE),meshSHA256=MESH,coldSHA256=COLD,
                    prepared=prep.relative_to(ROOT).as_posix(),preparedSHA256=sha(prep/"prepared.json"),
                    programs=native.relative_to(ROOT).as_posix(),steady=steady.relative_to(ROOT).as_posix(),
                    steadyFrame=base,steadyThermalSHA256=sha(steady/"thermal.json"),steadyBuild=sha(native/"steady/build.json"),
                    timeBuild=sha(native/"time/build.json"),transferSHA256=TRANSFER)
        phase(job,"physical","原生新算：从已通过收支核对的 14D 基流生成扰动，分别推进一个整步和两个半步")
        j=setup(work,native/"time",prepared=prep,prepared_sha=origin["preparedSHA256"],steady=steady/stages["stage"],
                base_frame=base,steady_build=native/"steady")
        assert j["changed"]==414 and abs(j["peakDeltaOverU"]-.01)<1e-13
        origin["seed"]=j["seed"];write_json(job/"construction.json",origin)
        run_steps(work,prepared=prep,plan=PLAN)
    else:
        work=completed.resolve();work.relative_to(ROOT/"outputs");origin=read(work.parent/"construction.json")
        receipt=read(work.parent/"completed.json")
        assert Path(receipt["source"]).resolve()==work and receipt["constructionSHA256"]==sha(work.parent/"construction.json")
        prep=directory(ROOT,origin["prepared"]);steady=directory(ROOT,origin["steady"]);base=origin["steadyFrame"]
        native=directory(ROOT,origin["programs"])
        assert origin["schema"]=="short65-cold-start-v1" and origin["meshSHA256"]==MESH==sha(prep/"mesh.su2")
        assert origin["caseSHA256"]==sha(CASE) and origin["preparedSHA256"]==sha(prep/"prepared.json")==PREP_KEY
        assert origin["transferSHA256"]==TRANSFER and origin["coldSHA256"]==sha(prep/"initial.csv")==COLD
        assert origin["steadyBuild"]==sha(native/"steady/build.json") and origin["timeBuild"]==sha(native/"time/build.json")
        assert read(work/"job.json")["baseFrame"]==base==check_steady(steady,prep,native)["frame"]
        assert origin["seed"]==sha(work/"seed.csv")==read(work/"job.json")["seed"] and origin["steadyThermalSHA256"]==sha(steady/"thermal.json")
        for branch,step in (("full",1),("half",1),("half",2)):
            prev=work/"seed.csv" if step==1 else work/"half/step1/restart_00001.csv"
            end=inspect_step(work,branch,step,prev,False,prepared=prep);inspect_step(work,branch,step,end,True,prepared=prep)
            endpoint(work,branch,step,prepared=prep)
        write_json(job/"construction.json",origin)
    phase(job,"mapping","逐面核对：输出五区壁压、热流、温度、近壁流向和七站对照")
    result=export(work,job/"export",prepared=prep,prepared_sha=origin["preparedSHA256"],base_frame=base,
                  steady_surface=steady.parent/"surface/surface.zip")
    assert result["faces"]==32480 and result["meshSHA256"]==MESH and result["steadyFrame"]==base and result["seed"]==origin["seed"]
    write_json(job/"completed.json",dict(pairSHA256=sha(job/"export/pair.json"),constructionSHA256=sha(job/"construction.json"),
               fullFrame=result["fullFrame"],halfFrame=result["halfFrame"],source=str(work),faces=result["faces"],
               areaM2=result["areaM2"],elapsedS=result["elapsedS"]))
    phase(job,"complete","原生三步、能量收支与五区逐面结果已保存")
