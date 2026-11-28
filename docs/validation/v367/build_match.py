"""复用已验证的两个观测对象，单独编译读入基态的本地回放入口。"""
import argparse
from pathlib import Path
import shutil
from build_observer import BASE,BIN,BUILD,DOC,invoke,replace,save,sha
import json


def build(src,out):
    src=src.resolve();out=out.resolve();src.relative_to(BASE);out.relative_to(BASE)
    old=json.loads((src/"build.json").read_bytes())
    assert sha(src/"tail_observe.exe")==old["executableSHA256"] and sha(src/"assembly.hpp")==old["observerSHA256"]
    out.mkdir(parents=True,exist_ok=False)
    for name in ("assembly.hpp","observed_upwind.hpp"):
        shutil.copyfile(src/name,out/name)
    shutil.copyfile(DOC/"match.hpp",out/"match.hpp")
    source=(src/"CNewtonIntegration.cpp").read_text(encoding="utf8")
    source=replace(source,'#include "assembly.hpp"','#include "assembly.hpp"\n#include "match.hpp"')
    marker="solvers[FLOW_SOL]->Set_OldSolution();"
    source=replace(source,marker,"TailBudget::bind_state(geometry,solvers[FLOW_SOL],config);\n  "+marker)
    unit=out/"CNewtonIntegration.cpp";unit.write_text(source,encoding="utf8")
    def command(name):
        r=json.loads((src/f"{name}_process.json").read_bytes())
        assert r["exitCode"]==0 and not r["livePids"] and Path(r["args"][0])==BIN/"g++.exe"
        return r["args"]
    args=command("CNewtonIntegration");obj=out/"CNewtonIntegration.obj"
    for key,val in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):args[args.index(key)+1]=str(val)
    args[-1]=str(unit)
    invoke("CNewtonIntegration",args,BUILD,out,180,3)
    args=command("link");args[args.index(str(src/"CNewtonIntegration.obj"))]=str(obj)
    exe=out/"tail_observe.exe";args[2]=str(exe)
    invoke("link",args,BUILD,out,180,3)
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<300*1024**2
    save(out/"build.json",dict(executableSHA256=sha(exe),baseBuildSHA256=sha(src/"build.json"),inputHookSHA256=sha(DOC/"match.hpp"),
         sourceSHA256=sha(Path(__file__)),newtonSourceSHA256=sha(unit),newtonObjectSHA256=sha(obj),
         reusedObjects={name:sha(src/name) for name in ("CIntegration.obj","CNumericsSIMD.obj")}))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("src",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();build(a.src,a.out)
