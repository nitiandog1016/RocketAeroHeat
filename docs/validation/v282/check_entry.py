"""独立重建生产入口的六套输入，与已经完成的原生运行逐项对照，不再推进流场。"""
from pathlib import Path
import json
import os
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
import beadSurfaceResult
import runThermistorCase
import thermistorBead
import thermistorBeadSector
import thermistorReaction
import thermistorTemperature


def check(out):
    source = ROOT/"outputs/v282_validation"
    native = out/"native"
    originals = {native/"axis":source/"outer01", native/"sector":source/"sector_outer01",
                 native/"strong":source/"strong_outer01"}
    for step in ("reaction0","plus20","reaction20"):
        originals[native/"probes"/step] = source/"probe_outer01"/step
    stages = {}
    # 保留真实的配置生成逻辑。每一阶段使用对应的已审查原生终场作为下阶段的输入。
    def compare_stage(folder):
        reference = originals[folder]
        files = ["case.cfg","mesh.su2"]
        files.extend(p.name for p in folder.iterdir()
                     if p.name in ("initial.npz","initial_flow.csv","nodes.npz","case.json","body.json","geometry.json"))
        pairs = {name:thermistorBead.sha(folder/name) for name in files}
        assert all(thermistorBead.sha(reference/name)==digest for name,digest in pairs.items())
        stages[folder.name] = pairs
    def use_final(prepare):
        return lambda parent,child: prepare(originals[parent],child)
    observe = thermistorReaction.observe
    def probe(parent,child,release):
        with patch.object(thermistorReaction,"terminal",side_effect=lambda f: terminal(originals.get(f,f))):
            observe(originals[parent],child,release=release)
    terminal = thermistorReaction.terminal
    collect = beadSurfaceResult.collect
    def read_final(parent,probes,result):
        return collect(originals[parent],originals[probes/"reaction0"].parent,result)
    with patch.object(thermistorBead,"solve",side_effect=compare_stage), \
         patch.object(thermistorBeadSector,"prepare",side_effect=use_final(thermistorBeadSector.prepare)), \
         patch.object(thermistorTemperature,"prepare",side_effect=use_final(thermistorTemperature.prepare)), \
         patch.object(thermistorReaction,"solve",side_effect=compare_stage), \
         patch.object(thermistorReaction,"observe",side_effect=probe), \
         patch.object(beadSurfaceResult,"collect",side_effect=read_final):
        runThermistorCase.run(native)
    assert len(stages)==6
    wall_hash=thermistorBead.sha(native/"result/wall.csv")
    assert wall_hash==thermistorBead.sha(source/"result_outer01/wall.csv")
    # 从真实入口检验关闭请求在启动引擎前终止任务。
    flag=Path(os.environ["RAH_THERMISTOR_CANCEL_FILE"]);flag.touch(exist_ok=False)
    try:
        thermistorBead.solve(native/"axis")
    except RuntimeError as exc:
        assert "用户已关闭" in str(exc)
    else:
        raise AssertionError("关闭后仍启动了求解器")
    assert not (native/"axis/engine_process.json").exists()
    report=dict(inputReplayOnly=True,nativeStepsRepeated=0,stages=stages,wallSHA256=wall_hash,
                cancellationBeforeSpawn=True)
    (out/"entry.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    check(Path(sys.argv[1]).resolve())
