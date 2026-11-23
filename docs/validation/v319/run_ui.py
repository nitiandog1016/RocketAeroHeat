"""从实际侧栏运行本轮验证，前端使用 v318 已发布的离线 WebGL 文件。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from runSU2NativeDebugger import observe, sha


def quoted(path):
    return "'" + path.as_posix().replace("'", "''") + "'"


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("name")
    p.add_argument("mode", choices=["fields", "native", "text"])
    a = p.parse_args()
    if not a.name.isidentifier():
        p.error("请指定新的运行目录名")
    out = ROOT / "outputs/v319_validation" / a.name
    out.mkdir(parents=True, exist_ok=False)
    front = out / "frontend"
    fixed = ["web", "createAerothermalWebTab.m", "createAerothermalGeometryView.m", "exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v318",*fixed],cwd=ROOT))) as z:
        z.extractall(front)
    files = [*front.rglob("*"), *HERE.glob("*.m"), HERE/"PLAN.md", Path(__file__), *ROOT.glob("*Thermistor*.m"),
             ROOT/"runAeroHeatSingleCase.m", ROOT/"rocketAeroHeatWorkbench.m", ROOT/"aeroHeatWorkflowState.m",
             ROOT/"buildAerothermalViewerData.m", ROOT/"exportWallPairEvidence.m", ROOT/"loadSTLMesh.m", ROOT/"readSTLTextCorners.m",
             ROOT/"outputs/v319_validation/fixture/probe-mm.stl", ROOT/"outputs/v319_validation/fixture/reference_ids.txt"]
    test=f"check_text_ui({quoted(out)})" if a.mode=="text" else f"check_import({quoted(out)},'{a.mode}')"
    command = f"cd({quoted(out)});addpath({quoted(ROOT)},{quoted(HERE)});addpath({quoted(front)},'-begin');{test}"
    args = [shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", command]
    if a.mode != "native":
        r = observe(args,out,{str(f):sha(f) for f in files if f.is_file()},360,300*1024**2)
        assert r["debuggerExitCode"]==0 and all(r["inputsStable"].values()) and not (r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"])
    else:
        # 原生子阶段自带独立的时间、内存和磁盘监督；外层保存工作台终态。
        import json
        import time
        locks = {str(f):sha(f) for f in files if f.is_file()}
        start = time.monotonic()
        with (out/"matlab.log").open("wb") as log:
            task = subprocess.Popen(args,cwd=str(out),stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
            (out/"invocation.json").write_text(json.dumps(dict(pid=task.pid,command=args,inputs=locks)),encoding="utf8")
            code = task.wait()
        stable = all(sha(Path(f))==h for f,h in locks.items())
        (out/"terminal.json").write_text(json.dumps(dict(pid=task.pid,exitCode=code,inputsStable=stable,seconds=time.monotonic()-start)),encoding="utf8")
        assert code==0 and stable
    assert (out/("text-ui.json" if a.mode=="text" else "ui.json")).is_file()
