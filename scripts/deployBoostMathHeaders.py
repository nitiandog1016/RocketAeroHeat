"""在E盘WSL安装固定签名Boost头文件，只服务成熟标量反算适配。"""
from pathlib import Path
import json,subprocess,time,shutil,os
from prepareOpenFoamControl import ROOT,save,sha


def main():
    from runOpenFoamProcess import run
    start=time.monotonic();out=ROOT/'outputs/v203_validation/boost01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    assert os.geteuid()==0
    package='libboost1.83-dev';version='1.83.0-2.1ubuntu3.2'
    meta=subprocess.check_output(['apt-cache','show',package+'='+version],text=True);(out/'package.txt').write_text(meta)
    save(out/'admission.json',dict(package=package,version=version,seconds=120,parentSeconds=150,engineBuilds=0,CFDRuns=0,
        workspace='E-drive WSL root filesystem',automaticContinuation=0))
    process=run(['apt-get','-y','--no-install-recommends','install',package+'='+version],out,out/'native',120,
        [Path(__file__),out/'package.txt',out/'admission.json',ROOT/'docs/validation/v203/BOUNDED_INVERSE_PLAN.md'],16,2048)
    actual=subprocess.check_output(['dpkg-query','-W','-f=${Version}',package],text=True)
    assert process['passed'] and actual==version
    sources=[]
    for name in ['math/tools/toms748_solve.hpp','math/tools/roots.hpp','version.hpp']:
        p=Path('/usr/include/boost')/name;q=out/'source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        sources.append(dict(path=str(p),sha256=sha(p)))
    license=Path('/usr/share/doc/libboost1.83-dev/copyright');shutil.copyfile(license,out/'copyright.txt')
    result=dict(passed=True,package=package,version=actual,process=process,sources=sources,licenseSHA256=sha(license),
        engineBuilds=0,CFDRuns=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(dict(passed=True,version=actual,elapsedSeconds=result['elapsedSeconds'])))


if __name__=='__main__':main()
