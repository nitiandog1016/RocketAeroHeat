"""独立核对v93失败前缀接续身份与实际计数，不重新求解或加载已清理矩阵。"""

from pathlib import Path
import hashlib
import json
import sys

import numpy as np
from scipy.io import loadmat
from scipy.sparse import issparse

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def identical(a, b):
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b)
        for key in a:
            identical(a[key], b[key])
    elif issparse(a):
        assert issparse(b) and a.shape == b.shape and (a != b).nnz == 0
    elif isinstance(a, (list, tuple)):
        assert isinstance(b, (list, tuple)) and len(a) == len(b)
        for left, right in zip(a, b):
            identical(left, right)
    elif isinstance(a, str):
        assert a == b
    else:
        np.testing.assert_array_equal(a, b)


def audit(folder):
    """模块一：直接比较保存数组、完整原行及来源SHA，不以布尔验收字段替代。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT)
    target = folder / "independent_recovery_identity.json"
    assert not target.exists()
    bindings = {}

    def record(path, expected=None):
        path = Path(path)
        digest = sha(path)
        assert expected is None or digest == expected, str(path)
        bindings[path.relative_to(ROOT).as_posix()] = digest
        return path

    def read(path):
        return json.loads(record(path).read_text(encoding="utf-8-sig"))

    recovery = read(folder / "recovery_binding.json")
    for key in ["failure", "context", "prefix", "preview", "independentProjectorProof"]:
        item = recovery[key]
        record(ROOT / item["file"], item["sha256"])
    failure = read(ROOT / recovery["failure"]["file"])
    assert failure["status"] == "terminal" and failure["exitCode"] == 1
    assert all(failure["sourceStableAfter"].values())
    process = read(folder / "coordinaterecover93_process.json")
    assert process["status"] == "terminal" and process["exitCode"] == 0
    assert all(process["sourceStableAfter"].values())
    step = folder / "continuous01/step01"
    prefix = loadmat(ROOT / recovery["prefix"]["file"], simplify_cells=True,
                     variable_names=["free", "freeCoordinates", "freeDetail", "originalRows"])
    active = loadmat(record(step / "active_direction.mat"), simplify_cells=True,
                     variable_names=["evidence", "directionReport", "preparationFile", "preparationSHA256"])
    context = loadmat(ROOT / recovery["context"]["file"], simplify_cells=True,
                      variable_names=["initialRows", "preparationFile", "preparationSHA256"])
    preview = loadmat(ROOT / recovery["preview"]["file"], simplify_cells=True,
                      variable_names=["value"])
    evidence, report = active["evidence"], active["directionReport"]
    for key in ["free", "freeCoordinates", "freeDetail"]:
        identical(prefix[key], evidence[key])
    identical(prefix["free"], preview["value"])
    identical(prefix["originalRows"], context["initialRows"])
    assert active["preparationFile"] == context["preparationFile"]
    assert active["preparationSHA256"] == context["preparationSHA256"]
    record(active["preparationFile"], active["preparationSHA256"])
    replay = report["freeCheckpointReplay"]
    assert replay["sha256"] == recovery["prefix"]["sha256"]
    assert replay["contextSHA256"] == recovery["context"]["sha256"]
    assert prefix["freeDetail"]["nativeIterations"] == replay["reusedNativeIterations"] == 128
    assert replay["newFreeNativeIterations"] == replay["newActualNSEvaluations"] == 0
    cycles = evidence["cycles"]
    cycles = [cycles] if isinstance(cycles, dict) else list(cycles)
    projected = sum(c["detail"]["nativeIterations"] for c in cycles)
    assert len(cycles) == 2 and projected == report["newNativeIterations"] == 256
    assert report["totalNativeIterations"] == projected + 128

    # 模块二：准备、探针与主线分开，接受一次并不代表稳态或物理解面积增加。
    row = read(step / "step.json")
    assert [row[k] for k in ["newPreparationActualNS", "newProbeActualNS", "newMainActualNS"]] == [1, 2, 4]
    assert row["newPreparationActualNSInRecovery"] == 0
    assert row["main"]["acceptedNonlinearIterations"] == 1
    record(step / "krylov01/wall_velocity_krylov_evidence.mat",
           "dc9fea4e454f137fa4e4ead736bed37ffb5407eed2cbc6b5f5012ca9e4ffd1ce")
    for rel in ["continuous01/independent_temperature_active.json",
                "independent_cell_interface_continuation.json", "row_real02/independent_row_conditioning.json"]:
        assert read(folder / rel)["status"].startswith("passed-")
    for name in ["solveRocketAugmentedQRDirection", "applyRocketAugmentedPatchQR", "applyRocketPackedGraph",
                 "applyRocketPackedGraphTranspose", "applyRocketPrimitiveUpdateMetric", "buildRocketAugmentedPatchQR"]:
        rel = "internalFlow/cutCell/" + name + ".m"
        record(ROOT / rel, failure["sourceSHA256"][rel])
    result = dict(status="passed-independent-same-state-exact-free-prefix-and-recovery-counter-identity",
                  inputSHA256=bindings, exactFreeAndCoordinatesEqual=True, exactOriginalRowsEqual=True,
                  reusedFreeIterations=128, newFreeIterations=0, projectedSolves=2, newProjectedIterations=256,
                  stageNSCounts=dict(preparation=1, probe=2, main=4), newNSCallsByThisAudit=0,
                  fullMatrixAuditReexecuted=False, prePruningNumericalAuditsPreservedAndHashed=True,
                  originalFailedExitRetained=1, recoveredExit=0, scriptSHA256=sha(__file__),
                  physicalSolution=False, newProductionSolvedAreaM2=0)
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(json.dumps({k: v for k, v in result.items() if k != "inputSHA256"}))


if __name__ == "__main__":
    audit(sys.argv[1])
