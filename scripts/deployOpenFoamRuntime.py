"""在E盘WSL环境部署固定官方OpenFOAM13；版本、签名来源和安装边界显式留档。"""
from pathlib import Path
import argparse,json,os,subprocess,time,urllib.request,hashlib

ROOT=Path('/mnt/e/RocketAeroHeat');BASE=ROOT/'outputs/v198_validation'
KEY=Path('/usr/share/keyrings/rah-openfoam.asc')
SOURCE=Path('/etc/apt/sources.list.d/rah-openfoam13.list')

def save(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def command(argv,out,seconds):
    """模块一：每项部署命令独立日志和明确时间上限，无shell拼接或全系统升级。"""
    start=time.monotonic()
    with out.with_suffix('.log').open('wb') as log:
        result=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=seconds,
            env={**os.environ,'DEBIAN_FRONTEND':'noninteractive','LC_ALL':'C'})
    row=dict(argv=argv,exitCode=result.returncode,elapsedSeconds=time.monotonic()-start)
    save(out.with_suffix('.json'),row);return row

def metadata():
    """模块二：仅增加具名专用仓库，固定官方签名密钥，先查实际包再决定安装。"""
    assert os.name=='posix' and os.geteuid()==0
    out=BASE/'openfoam_metadata01';out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    assert not KEY.exists() and not SOURCE.exists(),'不覆盖已有仓库配置'
    raw=urllib.request.urlopen('https://dl.openfoam.org/gpg.key',timeout=30).read()
    assert b'BEGIN PGP PUBLIC KEY BLOCK' in raw and len(raw)<100000
    (out/'openfoam.asc').write_bytes(raw);KEY.write_bytes(raw)
    line='deb [arch=amd64 signed-by=/usr/share/keyrings/rah-openfoam.asc] http://dl.openfoam.org/ubuntu noble main\n'
    SOURCE.write_text(line);(out/'source.list').write_text(line)
    commands=[]
    for name,argv in [('key',['gpg','--batch','--show-keys','--with-fingerprint',str(KEY)]),
        ('update',['apt-get','-o','Dir::Etc::sourcelist='+str(SOURCE),'-o','Dir::Etc::sourceparts=-',
            '-o','APT::Get::List-Cleanup=0','-o','Acquire::Retries=0','-o','Acquire::http::Timeout=25','update']),
        ('policy',['apt-cache','policy','openfoam13']),('package',['apt-cache','show','openfoam13'])]:
        r=command(argv,out/name,max(1,120-(time.monotonic()-start)));commands.append(r);assert r['exitCode']==0,name
    save(out/'report.json',dict(passed=True,keySHA256=hashlib.sha256(raw).hexdigest(),source=line,
        commands=commands,elapsedSeconds=time.monotonic()-start,CFDRuns=0,engineBuilds=0))
    print(json.dumps(dict(passed=True,elapsedSeconds=time.monotonic()-start)))

def install():
    """模块三：固定20260624发行包，先保存依赖方案再下载、核对包SHA和安装。"""
    assert os.name=='posix' and os.geteuid()==0
    out=BASE/'openfoam_install01';out.mkdir(exist_ok=False);start=time.monotonic()
    version='20260624';digest='75cb48bb7aee5907fb449ff4b520cbbef5e095c867d8c7320a6ca976aca62791'
    package='openfoam13='+version
    query=subprocess.check_output(['dpkg-query','-W','-f=${Package}\t${Version}\n'],text=True)
    (out/'packages_before.txt').write_text(query)
    base=['apt-get','--no-install-recommends','-o','Acquire::Retries=0','-o','Acquire::http::Timeout=30']
    simulation=command(base+['--assume-no','install',package],out/'dependency_plan',60)
    log=(out/'dependency_plan.log').read_text()
    assert simulation['exitCode']==1 and '151 MB of archives' in log and '745 MB of additional disk space' in log
    assert '0 to remove' in log and '1 upgraded, 44 newly installed' in log
    save(out/'admission.json',dict(version=version,debSHA256=digest,maximumSeconds=900,
        maximumDownloadBytes=1024**3,maximumInstalledBytes=4*1024**3,
        plannedDownloadMB=151,plannedInstalledMB=745,optionalParaView=False,CFDRuns=0))
    runs=[simulation]
    for name,argv in [('download',base+['-y','--download-only','install',package]),
                      ('install',base+['-y','--no-download','install',package])]:
        if name=='install':
            deb=Path('/var/cache/apt/archives/openfoam13_'+version+'_amd64.deb')
            assert deb.stat().st_size==118628390 and hashlib.sha256(deb.read_bytes()).hexdigest()==digest
        result=command(argv,out/name,max(1,900-(time.monotonic()-start)));runs.append(result)
        assert result['exitCode']==0,name
    installed=subprocess.check_output(['dpkg-query','-W','-f=${Version}','openfoam13'],text=True)
    assert installed==version
    (out/'packages_after.txt').write_text(subprocess.check_output(['dpkg-query','-W','-f=${Package}\t${Version}\n'],text=True))
    save(out/'report.json',dict(passed=True,version=version,debSHA256=digest,commands=runs,
        elapsedSeconds=time.monotonic()-start,CFDRuns=0,engineBuilds=0,rebootRequired=False))
    print(json.dumps(dict(passed=True,version=version,elapsedSeconds=time.monotonic()-start)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['metadata','install']);a=p.parse_args()
    metadata() if a.stage=='metadata' else install()
