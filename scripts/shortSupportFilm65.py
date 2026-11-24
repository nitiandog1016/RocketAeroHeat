"""由 Test 65 球面功率反求等温壁温，每次原生试算按候选温度独立计算膜温输运参数。"""
import argparse
import os
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

from shortSupportPower65 import CASE, PLAN as FIXED_PLAN, ROOT, TRANSFER, atomic_write_json, exclusive, read, sha, sources, trial, write_json
from shortSupportFilmProperties import film_state

PLAN=ROOT/"docs/validation/v341/PLAN.md"
MESH="a56a4c3dee8fd8bb74bc23ae4e85b71eaf7d5552dd99fccf8a6c8cb6a40789cd"
TOLERANCE=.0002


def inputs(job):
    """现有安装取历史流场的动量；没有历史资产时用试算温度从来流建立新场。"""
    warm=sources(job,PLAN,allow_bootstrap=False)
    if warm is not None:return warm
    import installGmshRuntime
    from installShort65Kernel import HOME_DIR, install, invoke
    from buildShort65Programs import build
    prep=job/"geometry";native=job/"programs"
    installGmshRuntime.main()
    invoke("mesh",[sys.executable,ROOT/"scripts/shortSupportMesh.py",prep,"--split-near","--split-sectors","--uniform-excess","20"],
           ROOT,job,180,2,bin_directory=Path(sys.executable).parent)
    assert sha(prep/"mesh.su2")==MESH
    if not (HOME_DIR/"kernel.json").is_file():install()
    build(native,prep)
    engine=native/"steady/shared.exe";release=native/"steady/wire.exe";initial=prep/"initial.csv"
    os.environ["PATH"]=str(HOME_DIR/"toolchain/w64devkit/bin")+os.pathsep+os.environ["PATH"]
    write_json(job/"source.json",dict(prepared=str(prep),preparedSHA256=sha(prep/"prepared.json"),frame=str(initial),
        frameSHA256=sha(initial),buildSHA256=sha(engine.parent/"build.json"),engineSHA256=sha(engine),releaseSHA256=sha(release),
        caseSHA256=sha(CASE),transferSHA256=sha(TRANSFER),planSHA256=sha(PLAN),uniformContractSHA256=sha(FIXED_PLAN.with_name("NATIVE.md"))))
    return prep,initial,engine,release


def solve(job):
    job=job.resolve();job.relative_to(ROOT/"outputs")
    if shutil.disk_usage(ROOT).free<16*1024**3:raise OSError("膜温闭合工况需要 16 GiB 可用磁盘。")
    job.mkdir(parents=True,exist_ok=False)
    with exclusive(job):
        base,frame,engine,release=inputs(job)
        case=read(CASE);target=case["observedBtuS"]*1055.05585262/case["supportAreaFactor"]
        history=[];records={};seed_excess=None
        def evaluate(name,excess):
            nonlocal frame,seed_excess
            state=trial(job,name,excess,(base,frame,engine,release),film=True,seed_excess=seed_excess,plan=PLAN)
            _,transport=film_state(case,excess)
            assert state["filmTransport"]==transport and state["temperatureK"]==transport["wallK"]
            history.append(name);records[name]=state
            frame=job/name/"run/early/initial_flow.csv";seed_excess=excess
            assert sha(frame)==state["earlyFrame"]
            atomic_write_json(job/"acceptedTrials.json",dict(order=history,trials=records))
            return state
        low=evaluate("low",20.);high=evaluate("high",40.)
        for iteration in range(1,4):
            qlo=low["earlyRegionW"][0];qhi=high["earlyRegionW"][0]
            if not qlo<target<qhi:raise ValueError("实算球面功率未夹住量热目标。")
            excess=low["excessK"]+(target-qlo)*(high["excessK"]-low["excessK"])/(qhi-qlo)
            assert low["excessK"]<excess<high["excessK"]
            name=f"root{iteration:02d}";state=evaluate(name,excess)
            if all(abs(state[k][0]/target-1)<TOLERANCE for k in ("earlyRegionW","lateRegionW")):break
            if state["earlyRegionW"][0]<target:low=state
            else:high=state
        else:raise RuntimeError("三个内部候选用尽，未闭合两端点球面功率。")
        out=job/name/"surface";s=read(out/"summary.json");archive=job/"power65.zip"
        with ZipFile(archive,"x") as z:
            for member in ("summary.json","wall.csv","stations.csv"):
                z.writestr(ZipInfo(member),(out/member).read_bytes(),compress_type=ZIP_DEFLATED,compresslevel=9)
        atomic_write_json(job/"completed.json",dict(schema="short65-sphere-power-film-v1",source=read(job/"source.json"),
            order=history,trials=records,balanceName=name,balance=state,targetSphereW=target,temperatureK=state["temperatureK"],
            faces=s["faces"],areaM2=s["areaM2"],experimentalJointAreaM2=0.,archiveSHA256=sha(archive)))
        atomic_write_json(job/"progress.json",dict(stage="complete",message="球面量热、候选壁温和原生输运参数均已闭合"))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("job",type=Path);solve(parser.parse_args().job)
