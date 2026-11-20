"""把已经通过第一半步验收的端点观察器绑定到第二个物理步和实际保存迭代。"""
import argparse
import gzip
import json
import os
from pathlib import Path

from run_step import HERE, ROOT, RUN, read, replace, sha
from observe_half import BUILD, invoke
from thermistorStartupSurface import write


def build(area):
    proof=read(area/"second_state/verified.json")
    mid=read(area/"first_energy/accepted.json")
    assert proof["passed"] and proof["step"]==2 and mid["accepted"]
    assert proof["oldFrameSHA256"]==mid["energy"]["acceptedFrameSHA256"]
    prior=RUN/"endpoint_build02";build_id=read(prior/"build.json")
    assert build_id["observerSHA256"]==mid["energy"]["observerSHA256"]
    source=prior/"CIntegration.cpp"
    assert sha(source)==build_id["sourceSHA256"]
    n=proof["innerRows"]
    old="config->GetTimeIter()==1 && config->GetInnerIter()==260 &&\n      config->GetnInner_Iter()==261"
    new=f"config->GetTimeIter()==2 && config->GetInnerIter()=={n} &&\n      config->GetnInner_Iter()=={n+1}"
    cpp=replace(source.read_text("utf8"),old,new)
    cached=json.loads(gzip.decompress((HERE.parent/"v272/evidence/build.json.gz").read_bytes()))["cachedBuildInputs"]
    assert all(sha(Path(p))==h for p,h in cached.items())
    out=area/"build01";out.mkdir(exist_ok=False)
    (out/source.name).write_text(cpp,encoding="utf8")
    for stage,limit in (("compile",120),("link",60)):
        recipe=read(prior/(stage+"_process.json"))
        assert recipe["exitCode"]==0 and Path(recipe["cwd"])==BUILD
        args=[str(out/Path(s).name) if s.startswith(str(prior)+os.sep) else s for s in recipe["args"]]
        assert sha(Path(args[0]))==build_id["compilerSHA256"]
        invoke(stage,args,BUILD,out,limit,3)
    assert all(sha(Path(p))==h for p,h in cached.items())
    assert sum(p.stat().st_size for p in out.iterdir())<300*1024**2
    write(out/"build.json",dict(observerSHA256=sha(out/"observe_bead.exe"),sourceSHA256=sha(out/source.name),
          firstHalfSourceSHA256=sha(source),compilerSHA256=build_id["compilerSHA256"],readoutInner=n,step=2))


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--root",type=Path,default=RUN/"pair")
    a=cli.parse_args();folder=a.root.resolve();folder.relative_to(RUN)
    build(folder)
