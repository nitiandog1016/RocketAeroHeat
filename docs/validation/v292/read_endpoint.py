"""补测已接受热敏珠端点的标量壁反力。附加装配计算方程，各状态量均不再更新。"""
import argparse
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[name]="1"
import numpy as np
from run_step import HERE,LIMIT,PREVIOUS,ROOT,RUN,input_file,observe,read,settings,sha
from observe_half import BUILD,invoke,sample
from check_half import RMS,columns,energy
from thermistorStartupSurface import checked,write


def build_probe(prior,out):
    original=prior/"CIntegration.cpp"
    assert sha(original)==read(prior/"build.json")["sourceSHA256"]
    text=original.read_text("utf8")
    start,rest=text.split("void CIntegration::Time_Integration(",1)
    marker="  switch (config->GetKind_TimeIntScheme()) {"
    assert rest.count(marker)==1
    guard=r'''
  // The appended evaluation reads the actual endpoint without advancing either system.
  if (config->GetTimeIter()==1 && config->GetInnerIter()==260 &&
      config->GetnInner_Iter()==261 && config->GetnMGLevels()==0 &&
      config->GetTime_Marching()==TIME_MARCHING::DT_STEPPING_1ST &&
      geometry->GetnPointDomain()==148274 && omp_get_num_threads()==1) {
    if (MainSolver!=FLOW_SOL && MainSolver!=SPECIES_SOL) std::abort();
    auto* log=std::fopen("readout_updates.csv","ab");
    if (!log) std::abort();
    std::fprintf(log,"%lu,%lu,%s\n",(unsigned long)config->GetTimeIter(),
                 (unsigned long)config->GetInnerIter(),MainSolver==FLOW_SOL?"flow":"scalar");
    std::fclose(log);
    return;
  }
'''
    (out/original.name).write_text(start+"void CIntegration::Time_Integration("+rest.replace(marker,guard+marker),encoding="utf8")
    import gzip,json
    cached=json.loads(gzip.decompress((HERE.parent/"v272/evidence/build.json.gz").read_bytes()))["cachedBuildInputs"]
    assert all(sha(Path(p))==h for p,h in cached.items())
    for stage,limit in (("compile",120),("link",60)):
        recipe=read(prior/(stage+"_process.json"))
        assert recipe["exitCode"]==0 and Path(recipe["cwd"])==BUILD
        args=[str(out/Path(s).name) if s.startswith(str(prior)+os.sep) else s for s in recipe["args"]]
        assert sha(Path(args[0]))==read(prior/"build.json")["compilerSHA256"]
        invoke(stage,args,BUILD,out,limit,3)
    assert all(sha(Path(p))==h for p,h in cached.items())
    assert sum(p.stat().st_size for p in out.iterdir())<300*1024**2
    write(out/"build.json",dict(observerSHA256=sha(out/"observe_bead.exe"),sourceSHA256=sha(out/original.name),
          originalObserverSourceSHA256=sha(original),compilerSHA256=read(prior/"build.json")["compilerSHA256"],
          cachedInputsVerified=len(cached)))


def audit(area,out):
    src=area/"first01";obs=area/"endpoint_observe";ctrl=area/"first_release"
    proof=area/"first_state/verified.json"
    assert read(proof)["prefixRows"]==240 and read(proof)["innerRows"]==260
    checked(obs)
    prior=(area/"first_observe/history_00001.csv").read_bytes().splitlines(keepends=True)
    current=(obs/"history_00001.csv").read_bytes().splitlines(keepends=True)
    assert len(prior)==len(current)==262 and prior[:261]==current[:261]
    stages=(obs/"readout_updates.csv").read_text("ascii").splitlines()
    assert stages==["1,260,flow","1,260,scalar"],stages
    data=columns(obs/"before_wall.csv",["id","p","u","v","w","theta"])
    data=data[np.argsort(data[:,0])]
    scratch=columns(obs/"restart_00001.csv",["PointID","Pressure","Velocity_x","Velocity_y","Velocity_z","Species_0"])
    np.testing.assert_array_equal(scratch[:,0],data[:,0])
    scratch_error=np.max(abs(scratch[:,1:]-data[:,1:]),axis=0)
    assert scratch_error.max()<2e-12,scratch_error
    energy(1,src,out,proof,obs,ctrl,settings,LIMIT)
    record=read(out/"energy.json")
    record.update(prefixRows=240,observerPrefixByteIdenticalRows=260,skippedAppendedUpdates=stages,
                  scratchFrameDifference=scratch_error.tolist(),acceptedPhysicalSteps=1)
    write(out/"accepted.json",record)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("build","observe","audit"))
    cli.add_argument("--root",type=Path,default=RUN);cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN)
    built=area/"endpoint_build02"
    out=built if a.stage=="build" else area/("endpoint_"+a.stage)
    if a.worker:
        assert a.stage=="audit"
        audit(area,out)
    else:
        out.mkdir(parents=True,exist_ok=False)
        common=[HERE/"PLAN.md",HERE/"ENDPOINT_READOUT.md",Path(__file__).resolve(),HERE/"run_step.py",
                PREVIOUS/"check_half.py",PREVIOUS/"observe_half.py",ROOT/"scripts/thermistorStartupSurface.py",
                area/"assembly_check/assembly.json",area/"first_state/verified.json"]
        if a.stage=="build":
            build_probe(area/"build01",out)
        elif a.stage=="observe":
            sample("observe",area/"first01",area/"first_state/verified.json",out,
                   built,settings,LIMIT,common)
        else:
            files=common+[area/"first_release/restart_00001.csv",area/"first01/restart_00001.csv",
                          area/"first01/seed_00000.csv",area/"first_observe/history_00001.csv"]
            files.extend(area/"endpoint_observe"/s for s in ("debugger_process.json","readout_updates.csv","history_00001.csv",
                "stages.csv","restart_00001.csv","before_wall.csv","after_wall.csv"))
            rec=observe([sys.executable,str(Path(__file__).resolve()),"audit","--root",str(area),"--worker"],out,
                        {str(p):sha(p) for p in files},120,100*1024**2)
            assert not rec["stopReason"] and rec["debuggerExitCode"]==0,out/"debugger.log"
