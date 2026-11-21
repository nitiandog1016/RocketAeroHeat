"""在已登记的带线热敏珠体网上，独立冷启动一次有限斜率的标量输运。"""
import argparse
import gzip
import hashlib
import json
import os
import sys
import uuid
import zipfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runThermistorLeadCase import configuration
from thermistorBead import ENGINE, ENGINE_SHA, read, sha, solve, write_json
from thermistorReaction import terminal

BASE=ROOT/"outputs/v300_validation/wired01"
KIT=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build"))
COMMON=("mesh.su2","body.json","case.json","nodes.npz","geometry.json","initial_flow.csv")


def candidate_config():
    # 限制器从冷态持续到延续和单步读数；流场仍用自己的原有重构配置。
    original=configuration()
    anchor="MUSCL_SPECIES= NO\n"
    assert original.count(anchor)==1
    return original.replace(anchor,"MUSCL_SPECIES= YES\nSLOPE_LIMITER_SPECIES= BARTH_JESPERSEN\n"
                            "LIMITER_ITER= 999999\nRAMP_MUSCL= NO\n")


def native_sources():
    lock=read(ROOT/"third_party/su2/source_lock.json")
    binary=ENGINE.parents[1]/"official-download.zip"
    assert sha(binary)==lock["archive"]["sha256"] and sha(ENGINE)==ENGINE_SHA
    with zipfile.ZipFile(binary) as z:
        member,=[n for n in z.namelist() if n.split("/")[-1]=="SU2_CFD.exe"]
        assert hashlib.sha256(z.read(member)).hexdigest()==ENGINE_SHA
    build=json.loads(gzip.decompress((HERE.parent/"v272/evidence/build.json.gz").read_bytes()))
    source_sha=sha(KIT/"su2.zip")
    assert source_sha==build["sourceArchiveSHA256"]
    files=("Common/src/CConfig.cpp","SU2_CFD/include/limiters/CLimiterDetails.hpp",
           "SU2_CFD/include/solvers/CScalarSolver.inl","SU2_CFD/src/solvers/CSpeciesSolver.cpp",
           "SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp")
    sources={}
    with zipfile.ZipFile(KIT/"su2.zip") as z:
        for rel in files:
            member,=[n for n in z.namelist() if n.endswith("/"+rel)]
            local=KIT/"source"/rel
            assert local.read_bytes()==z.read(member)
            sources[rel]=sha(local)
    return dict(version=lock["version"],officialArchiveSHA256=sha(binary),executableSHA256=ENGINE_SHA,
                sourceArchiveSHA256=source_sha,files=sources)


def prepare(folder):
    base=read(HERE.parent/"v300/evidence/inputs.json")["inputs"]
    assert all(sha(BASE/k)==v for k,v in base.items())
    terminal(BASE)
    assert sha(BASE/"restart.csv")==read(HERE.parent/"v300/evidence/summary.json")["sourceFrameSHA256"]
    assert (BASE/"case.cfg").read_text("ascii")==configuration()
    source=native_sources()
    folder.mkdir(parents=True,exist_ok=False)
    for name in COMMON:os.link(BASE/name,folder/name)
    cfg=candidate_config()
    (folder/"case.cfg").write_text(cfg,encoding="ascii")
    files=[HERE/"PLAN.md",Path(__file__).resolve(),ROOT/"scripts/runThermistorLeadCase.py",
           ROOT/"scripts/thermistorReaction.py",ROOT/"scripts/thermistorBead.py",HERE.parent/"v282/observe_bead.py",
           ROOT/"scripts/runEmbeddedSU2.py",*(folder/k for k in COMMON)]
    request=dict(workDirectory=str(folder),configFile="case.cfg",meshFile="mesh.su2",threads=1,maximumSeconds=720,
                 jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in files])
    write_json(folder/"engine_request.json",request)
    write_json(folder/"inputs.json",dict(native=source,baseFrame=sha(BASE/"restart.csv"),
               unchanged={name:sha(folder/name) for name in COMMON},configSHA256=sha(folder/"case.cfg"),
               requestSHA256=sha(folder/"engine_request.json"),sourceFiles={p:sha(ROOT/p) for p in request["sourceFiles"]}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("prepare","solve"));parser.add_argument("folder",type=Path)
    args=parser.parse_args();folder=args.folder.resolve();folder.relative_to(ROOT/"outputs")
    if args.action=="prepare":prepare(folder)
    else:
        spec=read(folder/"inputs.json")
        assert all(sha(folder/k)==v for k,v in spec["unchanged"].items())
        assert all(sha(ROOT/k)==v for k,v in spec["sourceFiles"].items())
        assert sha(folder/"case.cfg")==spec["configSHA256"] and (folder/"case.cfg").read_text("ascii")==candidate_config()
        assert sha(ENGINE)==spec["native"]["executableSHA256"]==ENGINE_SHA
        assert sha(folder/"engine_request.json")==spec["requestSHA256"]
        solve(folder)
        terminal(folder)
