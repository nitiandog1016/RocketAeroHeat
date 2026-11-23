"""为 B 珠编译原生逐节点壁值入口及同阶段观察器，复用固定 SU2 8.6.0 内核。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import re
import sys
from zipfile import ZipFile

from thermistorBead import ROOT, read, sha, write_json

sys.path.insert(0,str(ROOT/"docs/validation/v305"))
from observe_layer import INTEGRATION, TOOL, arguments, invoke, stage_source

MAIN = "SU2_CFD/src/SU2_CFD.cpp"
APPLY = r'''
  // 逐一匹配真实壁节点，将给定温度与强标量（释放测量时为零通量）送入原生接口。
  const auto markers = driver->GetMarkerIndices();
  const auto wall = markers.at("wall");
  const auto count = driver->GetNumberMarkerNodes(wall);
  std::ifstream boundary("wire_boundary.dat");
  std::map<unsigned long,std::pair<double,double>> values;
  for (unsigned long row=0; row<count; ++row) {
    unsigned long id;
    double temperature, scalar;
    if (!(boundary >> id >> temperature >> scalar) || !std::isfinite(temperature) ||
        !std::isfinite(scalar) || temperature < 1 || temperature > 1.5 ||
        scalar < 0 || scalar > 0.5 || !values.emplace(id,std::make_pair(temperature,scalar)).second)
      throw std::runtime_error("Invalid thermistor wall row");
  }
  std::string extra;
  if (boundary >> extra) throw std::runtime_error("Extra thermistor wall row");
  for (unsigned long vertex=0; vertex<count; ++vertex) {
    const auto id=driver->GetNodeGlobalIndex(driver->GetMarkerNode(wall,vertex));
    const auto entry=values.at(id);
    driver->SetMarkerCustomTemperature(wall,vertex,entry.first);
    driver->SetMarkerCustomScalar(wall,vertex,{entry.second});
    values.erase(id);
  }
  std::cout << "Applied thermistor wall nodes: " << count << std::endl;
'''


def build(out):
    kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",TOOL)).resolve()
    source=kit/"source";work=kit/"build-v117";bin_dir=kit/"toolchain/w64devkit/bin"
    cache=json.loads(gzip.decompress((ROOT/"docs/validation/v272/evidence/build.json.gz").read_bytes()))
    old,=[Path(p) for p in cache["cachedBuildInputs"] if Path(p).as_posix().endswith("/source/"+INTEGRATION)]
    old_kit=old.parents[len(Path(INTEGRATION).parts)]
    locked={kit/Path(k).relative_to(old_kit):h for k,h in cache["cachedBuildInputs"].items()}
    assert all(sha(p)==h for p,h in locked.items())
    assert sha(kit/"su2.zip")==cache["sourceArchiveSHA256"] and sha(bin_dir/"g++.exe")==cache["compilerSHA256"]
    paths=[MAIN,INTEGRATION,"SU2_CFD/src/solvers/CSpeciesSolver.cpp","SU2_CFD/src/solvers/CIncNSSolver.cpp",
           "SU2_CFD/include/drivers/CDriverBase.hpp","SU2_CFD/src/drivers/CDriverBase.cpp"]
    with ZipFile(kit/"su2.zip") as archive:
        for rel in paths:
            member,=[n for n in archive.namelist() if n.endswith("/"+rel)]
            assert archive.read(member)==(source/rel).read_bytes()
    out.mkdir(parents=True,exist_ok=False)
    main=(source/MAIN).read_text("utf8")
    needle="  driver->StartSolver();"
    assert main.count(needle)==1
    entries=read(work/"compile_commands.json")
    outputs=[]
    for rel,code in ((MAIN,"#include <fstream>\n"+main.replace(needle,APPLY+needle)),
                     (INTEGRATION,stage_source((source/INTEGRATION).read_text("utf8")))):
        path=source/rel;copy=out/path.name
        copy.write_text(code,encoding="utf8")
        entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
        args=arguments(entry["command"]);args[0]=bin_dir/"g++.exe";obj=out/(path.stem+".obj")
        for flag,val in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):args[args.index(flag)+1]=val
        args[-1:]=["-iquote",path.parent,copy]
        invoke(path.stem,args,work,out,180,3,bin_directory=bin_dir);outputs.append(obj)
    link=re.search(r"^build SU2_CFD/src/SU2_CFD\.exe: cpp_LINKER ([^ |]+)[^\n]*\n LINK_ARGS = ([^\n]+)",(work/"build.ninja").read_text("utf8"),re.M)
    programs={}
    for count,name in ((1,"wire"),(2,"wire_observer")):
        exe=out/(name+".exe")
        invoke(name,[bin_dir/"g++.exe","-o",exe,*outputs[:count],*arguments(link[2])],work,out,180,3,bin_directory=bin_dir)
        programs[name]=sha(exe)
    assert all(sha(p)==h for p,h in locked.items())
    write_json(out/"build.json",dict(programs=programs,driverSourceSHA256=sha(Path(__file__)),releaseSHA256=cache["sourceArchiveSHA256"],
               compilerSHA256=cache["compilerSHA256"],sourceSHA256={rel:sha(source/rel) for rel in paths},
               generatedSHA256={p.name:sha(p) for p in out.glob("*.cpp")},cachedInputs=len(locked),license="LGPL-2.1-or-later"))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("out",type=Path)
    dest=parser.parse_args().out.resolve();dest.relative_to(ROOT/"outputs");build(dest)
