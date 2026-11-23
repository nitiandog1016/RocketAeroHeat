"""打包复核所需的小原始记录；全节点和重启场仅登记本地位置与 SHA。"""
import hashlib
import io
import json
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from verify_flux import BASE, DATA, HERE, OLD, ROOT, read, sha, sorted_by_id, table


def make_evidence(target):
    result = read(DATA/"supportedFlux65.json")
    detail = BASE/"detail01"
    small = []
    large = []
    for folder in (OLD/"run01/steady01", BASE/"run01/steady01", BASE/"scout01", detail):
        small += [folder/n for n in ("case.cfg", "debugger_process.json", "history.csv", "debugger.log")]
        large += list(folder.glob("restart_*.csv"))
        if (folder/"replay.json").exists():
            small += [folder/n for n in ("replay.json", "probe.in", "extrema.csv")]
            large.append(folder/"first_update.csv")
        else:
            small += [folder.parent/n for n in ("failure.json", "qualified.json")]
    small += [detail/n for n in ("flow_edges.csv", "scalar_edges.csv", "flow_far.csv", "scalar_far.csv", "row.csv", "matrix.csv")]
    small += list((BASE/"build01").glob("*.json"))
    small += [BASE/"ui01"/n for n in ("ui.json", "debugger_process.json")]
    small.append(OLD/"prep02/prepared/prepared.json")
    large += [OLD/"prep02/prepared"/n for n in ("mesh.su2", "initial.csv", "profile.dat", "mapping.npz", "transfer.npz")]
    large += [BASE/"build01/probe.exe", ROOT/"outputs/v317_validation/build01/wire.exe"]
    # 六个越界节点及目标节点矩阵邻点，保留足以离线重算局部残差行的坐标和原更新。
    states = sorted_by_id(table(detail, "first_update.csv", header=True))
    ids = np.union1d(result["first"]["belowIDs"], table(detail, "matrix.csv")[:, 0].astype(int))
    with np.load(OLD/"prep02/prepared/mapping.npz") as mapping:
        local = np.column_stack((states[ids, 0], mapping["xyz"][ids], states[ids, 1:]))
    subset = io.StringIO()
    np.savetxt(subset, local, delimiter=",", fmt=["%d"]+["%.17g"]*8, comments="",
               header="id,xD,yD,zD,theta_before,theta_after,spatial,delta,relax")
    def key(p):
        return p.relative_to(ROOT).as_posix()
    def digest(b):
        return hashlib.sha256(b).hexdigest()
    payload = {key(p): p.read_bytes() for p in small}
    payload["first_nodes.csv"] = subset.getvalue().encode("ascii")
    manifest = {"included": {k: digest(v) for k, v in payload.items()}, "largeLocal": {key(p): sha(p) for p in large}}
    payload["manifest.json"] = (json.dumps(manifest, indent=2)+"\n").encode("utf8")
    with ZipFile(target, "x", compression=ZIP_DEFLATED) as archive:
        for name, value in payload.items():
            archive.writestr(name, value)
    with ZipFile(target) as archive:
        assert set(archive.namelist()) == set(payload)
        assert all(archive.read(k) == v for k, v in payload.items())


if __name__ == "__main__":
    evidence = HERE/"evidence"
    evidence.mkdir(exist_ok=True)
    make_evidence(evidence/"receipts.zip")
