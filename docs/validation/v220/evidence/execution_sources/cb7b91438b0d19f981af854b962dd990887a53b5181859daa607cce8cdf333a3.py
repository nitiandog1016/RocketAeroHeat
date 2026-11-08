"""模块职责：生成实际依赖图和隔离薄归档；所有输出仅在v220目录。"""
from pathlib import Path
import json,subprocess,sys
from prepareSU2PhaseAirV220 import ROOT,BASE,BUILD,SOURCE,OLD_SOURCE,OLD_BUILD,CPP,read,save,sha

PREP=BASE/'prepare03'
def compile_commands():
    commands=read(PREP/'commands.json');result=[]
    for name in CPP:
        cmd=commands[name];obj=Path(cmd[cmd.index('-o')+1])
        result.append(cmd[:-2]+['-MD','-MF',str(obj)+'.d']+cmd[-2:])
    return result

def dependencies():
    folder=BASE/'dependencies01';folder.mkdir(exist_ok=False);paths={};commands=[]
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    changes=read(PREP/'changes.json')
    for i,basecmd in enumerate(compile_commands()):
        cmd=[];j=0
        while j<len(basecmd):
            if basecmd[j] in ('-o','-MF'):j+=2;continue
            if basecmd[j] in ('-c','-MD'):j+=1;continue
            cmd.append(basecmd[j]);j+=1
        dep=folder/('unit%d.d'%i);cmd+=['-M','-MT','dep','-MF',str(dep)];commands.append(cmd)
        subprocess.check_call(cmd,cwd=str(OLD_BUILD),timeout=20)
        text=dep.read_text(encoding='utf8').replace('\\\n',' ')
        assert text.startswith('dep:');names=text[4:].split()
        for name in names:
            p=Path(name)
            if not p.is_absolute():p=OLD_BUILD/p
            p=p.resolve();assert p.is_file(),str(p)
            digest=sha(p);paths[str(p)]=digest
            # 对原SU2源/新增隔离源分别绑定；外部编译器/Boost实际依赖完整锁定。
            try:relative=p.relative_to(SOURCE).as_posix()
            except ValueError:
                try:relative=p.relative_to(OLD_SOURCE).as_posix()
                except ValueError:continue
                assert tree[str(OLD_SOURCE/relative)]==digest,relative
            else:
                if relative in changes:assert changes[relative]['patchedSHA256']==digest
                elif str(OLD_SOURCE/relative) in tree:assert tree[str(OLD_SOURCE/relative)]==digest,relative
    save(folder/'commands.json',commands);save(folder/'inputs.json',paths)
    save(folder/'report.json',dict(passed=True,translationUnits=len(commands),inputFiles=len(paths),newBuild=False))
    print(json.dumps(dict(passed=True,dependencyFiles=len(paths))))

def link():
    from prepareSU2PhaseAirV220 import BOOST
    folder=BASE/'link01';folder.mkdir(exist_ok=False)
    selected=read(PREP/'selected_objects.json');old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for name,h in old.items():assert sha(OLD_BUILD/name)==h,name
    compiler=Path(compile_commands()[0][0]);ar=compiler.parent/'ar.exe';archives={};members_report={}
    for group in ('Common','SU2_CFD'):
        prior=OLD_BUILD/(group+'/src/'+('libSU2Common.a' if group=='Common' else 'libSU2core.a'))
        raw=subprocess.check_output([str(ar),'t',str(prior)],cwd=str(OLD_BUILD)).decode().splitlines()
        members=[]
        for n in raw:
            p=Path(n)
            if not p.is_absolute():p=OLD_BUILD/p
            name=p.resolve().relative_to(OLD_BUILD).as_posix()
            member=Path(selected.get(name,str(p)));assert member.is_file(),str(member);members.append(str(member))
        target=BUILD/prior.name
        assert not target.exists()
        subprocess.check_call([str(ar),'crsT',str(target)]+members,cwd=str(OLD_BUILD),timeout=10)
        archives[prior.name]=str(target);members_report[prior.name]=[dict(path=p,sha256=sha(Path(p))) for p in members]
    command=read(ROOT/'outputs/v167_validation/build01/commands.json')['link'];cmd=[]
    for i,x in enumerate(command):
        if i==2:cmd.append(str(BUILD/'SU2_CFD.exe'))
        elif x.endswith('.a'):cmd.append(archives[Path(x).name])
        elif x.endswith('.obj'):cmd.append(str(OLD_BUILD/x))
        else:cmd.append(x)
    save(folder/'members.json',members_report);save(folder/'command.json',cmd)
    subprocess.check_call(cmd,cwd=str(OLD_BUILD),timeout=25)
    save(folder/'report.json',dict(passed=True,executableSHA256=sha(BUILD/'SU2_CFD.exe'),overriddenObjects=len(selected),
        oldObjectsStable=True,sourceModel='rah-su2-phase-air-saneg-eplus300k-v1',nativeVerification=False,physicalQualification=False))

if __name__=='__main__':
    if sys.argv[1]=='dependencies':dependencies()
    elif sys.argv[1]=='link':link()
    else:raise ValueError('No CFD entry')
