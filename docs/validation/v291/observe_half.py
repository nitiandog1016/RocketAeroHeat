"""用固定的单对象构建配方，为热敏珠两个半步读取同时间层的强标量壁反力。"""
import argparse
import gzip
import json
import os
from pathlib import Path

from run_half import (BIN, ENGINE_SHA, HERE, MAX_INNER, OFFICIAL, ROOT, RUN, config,
                      input_file, observe, read, replace, sha)
from observe_heat import BUILD, FILES, SRC, instrument, invoke, save
from thermistorStartupSurface import checked


def build(out):
    reference=HERE.parent/"v290/evidence"
    base=read(reference/"build.json")
    archive=json.loads(gzip.decompress((reference/"processes.json.gz").read_bytes()))
    cached=json.loads(gzip.decompress((HERE.parent/"v272/evidence/build.json.gz").read_bytes()))["cachedBuildInputs"]
    assert all(sha(Path(p))==h for p,h in cached.items())
    original=SRC/FILES[0]
    assert sha(original)==base["originalSourceSHA256"]
    text=instrument(original.read_text("utf8"))
    text=replace(text,"const auto iteration=config->GetInnerIter();",
                 "const auto iteration=config->GetInnerIter();\n  const auto last=config->GetnInner_Iter()-1;")
    text=replace(text,"config->GetTimeIter()==1 && (iteration==199 || iteration==200)",
                 "(config->GetTimeIter()==1 || config->GetTimeIter()==2) && (iteration+1==last || iteration==last)")
    assert text.count("iteration==200 && dump")==2
    (out/original.name).write_text(text.replace("iteration==200 && dump","iteration==last && dump"),encoding="utf8")
    prior=str(ROOT/"outputs/v290_validation/build01")+os.sep
    for stage,limit in (("compile",120),("link",60)):
        recipe=archive["build01/"+stage+"_process.json"]
        assert Path(recipe["cwd"])==BUILD and recipe["exitCode"]==0
        args=[str(out/Path(s).name) if s.startswith(prior) else s for s in recipe["args"]]
        assert sha(Path(args[0]))==base["compilerSHA256"]
        invoke(stage,args,BUILD,out,limit,3)
    assert all(sha(Path(p))==h for p,h in cached.items())
    assert sum(p.stat().st_size for p in out.iterdir())<300*1024**2
    save(out/"build.json",dict(observerSHA256=sha(out/"observe_bead.exe"),sourceSHA256=sha(out/original.name),
                               compilerSHA256=base["compilerSHA256"],cachedInputsVerified=len(cached)))


def observer_config(index,rows):
    text=replace(config(index),f"INNER_ITER= {MAX_INNER}\n",f"INNER_ITER= {rows+1}\n")
    if rows<MAX_INNER:
        text=replace(text,"CONV_STARTITER= 50\n",f"CONV_STARTITER= {rows+2}\n")
    return text


def release_config(index):
    text=replace(config(index),f"INNER_ITER= {MAX_INNER}\n","INNER_ITER= 1\n")
    return replace(text,"(wall, VALUE, 1.0)","(wall, FLUX, 0.0)")


def sample(stage,source,proof_file,out,built):
    proof=read(proof_file)
    index=proof["step"]
    assert proof["passed"] and sha(source/f"restart_{index:05}.csv")==proof["frameSHA256"]
    checked(source)
    assert (source/"case.cfg").read_text("ascii")==config(index)
    old_name=f"seed_{index-1:05}.csv"
    locks=[HERE/"PLAN.md",HERE/"run_half.py",Path(__file__).resolve(),proof_file,
           source/"case.cfg",source/"debugger_process.json",ROOT/"scripts/runSU2NativeDebugger.py"]
    if stage=="observe":
        exe=built/"observe_bead.exe"
        assert sha(exe)==read(built/"build.json")["observerSHA256"]
        text=observer_config(index,proof["innerRows"])
        initial=source/old_name
        locks.extend([built/"build.json",built/"CIntegration.cpp"])
        limit,mb=420,350
    else:
        exe=OFFICIAL
        assert sha(exe)==ENGINE_SHA
        text=release_config(index)
        initial=source/f"restart_{index:05}.csv"
        limit,mb=90,180
    os.link(initial,out/old_name)
    os.link(input_file("mesh.su2"),out/"mesh.su2")
    (out/"case.cfg").write_text(text,encoding="ascii")
    locks.extend([exe,out/"case.cfg",out/"mesh.su2",out/old_name])
    os.environ["PATH"]=str(BIN)+os.pathsep+os.environ["PATH"]
    observe([str(exe),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in locks},limit,mb*1024**2)
    checked(out)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=["build","observe","release"])
    cli.add_argument("output",type=Path)
    cli.add_argument("--source",type=Path)
    cli.add_argument("--proof",type=Path)
    cli.add_argument("--build",type=Path)
    a=cli.parse_args()
    target=a.output.resolve()
    target.relative_to(RUN)
    target.mkdir(parents=True,exist_ok=False)
    if a.stage=="build":
        build(target)
    else:
        sample(a.stage,a.source.resolve(),a.proof.resolve(),target,a.build.resolve() if a.build else None)
