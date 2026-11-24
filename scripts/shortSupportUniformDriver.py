"""为 Short 65 的单一等温边界编译校验，复用原生标量方程、流动对象和固定库。"""
from pathlib import Path
import shutil

from build_shared import invoke, once, read, sha, write_json


def build(out,original,theta):
    if not 20/90<=theta<=40/90:
        raise ValueError("本量热试算的登记超温为 20–40 K。")
    out.mkdir()
    manifest=read(original/"build.json")
    for name in ("CSpeciesSolver.cpp","capacity.hpp","silver.hpp","far_mass.hpp","CScalarSolver.inl"):
        assert sha(original/name)==manifest["files"][name]
        if name!="CScalarSolver.inl":shutil.copyfile(original/name,out/name)
    key="    if(!wallCount || std::abs(upper-.4683950617283952)>1e-13) std::abort();"
    guard=f"""    constexpr double prescribedTheta = {theta:.17g};
    if (wallCount != 11906 || std::abs(upper-prescribedTheta)>1e-13) std::abort();
    for (const auto& node : nodesData) {{
      if (node.wall && (!std::isfinite(node.old) || std::abs(node.old-prescribedTheta)>1e-13)) std::abort();
    }}"""
    text=(original/"CScalarSolver.inl").read_text("utf8")
    (out/"CScalarSolver.inl").write_text(once(text,key,guard),encoding="utf8")
    compile_job=read(original/"CSpeciesSolver_process.json")
    link_job=read(original/"shared_process.json")
    for job in (compile_job,link_job):
        assert job["exitCode"]==0 and not(job["stopReason"] or job["livePids"])
    compile_args=compile_job["args"][:];cc=Path(compile_args[0])
    assert sha(cc)=="315c23bc298aba9eaf2aece174fa7fcb0e404247176663e6a4c7707a60f8e444"
    assert Path(compile_args[-1])==original/"CSpeciesSolver.cpp" and Path(link_job["args"][0])==cc
    obj=out/"CSpeciesSolver.obj"
    for flag,dest in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):
        compile_args[compile_args.index(flag)+1]=str(dest)
    compile_args[-1]=str(out/"CSpeciesSolver.cpp")
    cwd=Path(compile_job["cwd"]);assert Path(link_job["cwd"])==cwd
    invoke("species",compile_args,cwd,out,180,3,bin_directory=cc.parent)
    args=link_job["args"][:];exe=out/"uniform.exe"
    assert args.count(str(original/"CSpeciesSolver.obj"))==1 and Path(args[2])==original/"shared.exe"
    args[2]=str(exe);args[args.index(str(original/"CSpeciesSolver.obj"))]=str(obj)
    objects={str(Path(x)):sha(Path(x)) for x in args[3:] if x.endswith(".obj")}
    for file,digest in objects.items():
        path=Path(file)
        if path==obj:continue
        expected=manifest["files"].get(path.name,manifest.get("main") if path.name=="SU2_CFD.obj" else None)
        assert digest==expected
    invoke("link",args,cwd,out,180,3,bin_directory=cc.parent)
    write_json(out/"build.json",dict(expectedTheta=theta,wallNodes=11906,base=sha(original/"build.json"),
        template=sha(original/"CScalarSolver.inl"),uniformSource=sha(out/"CScalarSolver.inl"),objects=objects,exe=sha(exe)))
    return exe
