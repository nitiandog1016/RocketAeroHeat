"""用两版已保存的原生终态复查 7.5° 冷启动命令的每个输入及最终表。"""
import json
from pathlib import Path
import runpy
import sys
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
import beadSurfaceResult
import thermistorBead as bead
import thermistorBeadSector as sector
import thermistorReaction as reaction
import thermistorTemperature as temperature


def check(output):
    old=ROOT/"outputs/v282_validation";new=ROOT/"outputs/v283_validation"
    frozen={output/"axis":old/"outer01",output/"sector":new/"sector01",output/"strong":new/"strong01"}
    for name in ("reaction0","plus20","reaction20"):
        frozen[output/"probes"/name]=new/"probes01"/name
    matched={}
    def compare_inputs(folder):
        files=[p for p in folder.iterdir() if p.suffix in (".cfg",".su2",".npz",".csv",".json")
               and p.name!="engine_request.json"]
        assert len(files)>=3
        digests={p.name:bead.sha(p) for p in files}
        assert all(bead.sha(frozen[folder]/n)==h for n,h in digests.items())
        matched[folder.name]=digests
    prepare_sector,prepare_temperature=sector.prepare,temperature.prepare
    observe,terminal,collect=reaction.observe,reaction.terminal,beadSurfaceResult.collect
    def observe_saved(parent,folder,release):
        with patch.object(reaction,"terminal",lambda p:terminal(frozen.get(p,p))):
            return observe(frozen[parent],folder,release=release)
    # 运行公开命令的参数解析与阶段调度，只在实际引擎调用处对照已有输入字节。
    with patch.object(sys,"argv",["runThermistorCase.py",str(output),"--angular-cells","4"]), \
         patch.object(bead,"solve",compare_inputs), patch.object(reaction,"solve",compare_inputs), \
         patch.object(sector,"prepare",lambda p,f,n:prepare_sector(frozen[p],f,n)), \
         patch.object(temperature,"prepare",lambda p,f:prepare_temperature(frozen[p],f)), \
         patch.object(reaction,"observe",observe_saved), \
         patch.object(beadSurfaceResult,"collect",lambda p,g,f:collect(frozen[p],frozen[g/"reaction0"].parent,f)):
        runpy.run_module("runThermistorCase",run_name="__main__")
    assert len(matched)==6
    digest=bead.sha(output/"result/wall.csv")
    assert digest==bead.sha(new/"result01/wall.csv")
    (output/"entry.json").write_text(json.dumps(dict(nativeStepsRepeated=0,angularCells=4,matched=matched,
                                                    wallSHA256=digest),indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    check(Path(sys.argv[1]).resolve())
