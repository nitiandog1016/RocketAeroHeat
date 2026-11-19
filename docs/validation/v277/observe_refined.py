"""复用已审核的低速半域观测，只在本次 29,425 节点网格上回读同阶段能量账。"""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"docs/validation/v273"))
from observe_sector import build, sample
from check_sector_flux import audit
from cm011Tunnel import write_json
from readCM011Inclined import SCHEMA, state_at_frame


def main(action, base):
    base=base.resolve()
    base.relative_to(ROOT/"outputs")
    inputs=json.loads((base/"gui01/inputs.json").read_bytes())
    assert inputs["refinedNodes"]==29425
    if action=="build":
        build(base/"build01",domain_nodes=inputs["refinedNodes"])
        return
    run=Path((base/"gui01/run_folder.txt").read_text("utf-8-sig").strip())/"flow"
    assert json.loads((run/"prepared.json").read_bytes())==inputs["inputs"]
    if action=="sample":
        sample(run,base/"official01",None,state_at_frame)
        sample(run,base/"trace01",base/"build01",state_at_frame)
    elif action=="audit":
        result=audit(run,base/"official01",base/"trace01",state_at_frame,meridian_nodes=49,
                     copies=2,symmetry=("mirror",),schema=SCHEMA)
        write_json(base/"energy.json",result)
        assert result["belowHalfPercent"],"边界能量差额超出预定门槛"


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("build","sample","audit"))
    parser.add_argument("base",type=Path)
    args=parser.parse_args()
    main(args.action,args.base)
