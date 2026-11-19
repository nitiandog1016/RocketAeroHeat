"""复用三维只读装配观测，在 CM-011 的单镜面半域核对热收支。"""
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"docs/validation/v273"))
from observe_sector import build, sample
from check_sector_flux import audit
from cm011Tunnel import write_json
from readCM011Inclined import SCHEMA, state_at_frame


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("build","official","trace","audit"))
    parser.add_argument("run",type=Path)
    parser.add_argument("new_output",type=Path)
    parser.add_argument("--build",type=Path)
    parser.add_argument("--official",type=Path)
    parser.add_argument("--trace",type=Path)
    args=parser.parse_args()
    parent=args.run.resolve()
    out=args.new_output.resolve()
    parent.relative_to(ROOT/"outputs")
    out.relative_to(ROOT/"outputs")
    if args.action=="build":
        build(out,domain_nodes=15325)
    elif args.action=="audit":
        result=audit(parent,args.official.resolve(),args.trace.resolve(),state_at_frame,
                     meridian_nodes=49,copies=2,symmetry=("mirror",),schema=SCHEMA)
        write_json(out,result)
    else:
        sample(parent,out,args.build.resolve() if args.action=="trace" else None,state_at_frame)
