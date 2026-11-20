"""检查相对路径预检、网格错配拒绝和两档已确认体网的三维初态，不增加求解步。"""
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
import runThermistorCase
from thermistorBead import read,sha
from thermistorBeadSector import prepare


class ColdInputReached(Exception):
    """调用方已走过实际的预检和冷初场生成。"""


def check(out,workflow):
    out.mkdir(exist_ok=False)
    results={}
    for cells,key,earlier in [(2,"grid15","v282_validation/sector_outer01"),(4,"grid7p5","v283_validation/sector01")]:
        grid=Path(workflow[key]);digest=sha(grid/"prepared.json")
        for incorrect_n,incorrect_digest in [(6-cells,digest),(cells,"0"*64)]:
            rejected=out/f"reject_{cells}"
            try:
                runThermistorCase.run(rejected,incorrect_n,grid,incorrect_digest)
            except ValueError:
                assert not rejected.exists()
            else:
                raise AssertionError("wrong grid accepted")
        run=out/f"cold_{cells}"
        with patch.object(runThermistorCase.thermistorBead,"solve",side_effect=ColdInputReached):
            try:
                runThermistorCase.run(run,cells,Path(os.path.relpath(grid)),digest)
            except ColdInputReached:
                assert sha(run/"axis/mesh.su2")==sha(ROOT/"outputs/v282_validation/outer01/mesh.su2")
            else:
                raise AssertionError("cold stage missing")
        sector=out/f"sector_{cells}"
        prepare(ROOT/"outputs/v282_validation/outer01",sector,cells,grid,digest)
        reference=ROOT/"outputs"/earlier
        names=["mesh.su2","nodes.npz","initial_flow.csv","case.cfg","case.json","body.json","geometry.json"]
        assert all(sha(sector/n)==sha(reference/n) for n in names)
        results[str(cells)]={n:sha(sector/n) for n in names}
    (out/"check.json").write_text(json.dumps(dict(stages=results,refusedWrongResolutionAndDigest=True,
        relativePathValidated=True,nativeStepsRepeated=0),indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    check(Path(sys.argv[1]).resolve(),read(Path(sys.argv[2])))
