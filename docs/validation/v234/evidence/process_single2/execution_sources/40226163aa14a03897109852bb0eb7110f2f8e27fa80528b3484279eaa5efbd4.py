"""模块职责：复用原小壳与归档，登记并监督单/双黏性域整体控制；不运行真实CFD。"""
from pathlib import Path
import json,shutil,subprocess,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v234_validation';DOC=ROOT/'docs/validation/v234';BUILD=BASE/'build01'
CPP=ROOT/'native/su2PhaseAwareAir/multiZoneCompatibilityV234.cpp'
OLD=ROOT/'outputs/v221_validation';GAS=ROOT/'outputs/v220_validation'
OLD_SOURCE=ROOT/'tools/runtime/su2-thermal-build/source';OLD_BUILD=ROOT/'tools/runtime/su2-thermal-build/build-v117'
CLOSED=False


def split_mesh(out):
    """模块一：只分割现有40个母单元，保留原坐标、原单元内顺序和全部物理边界。"""
    p=OLD/'control01/mesh.su2';rows=p.read_text(encoding='utf8').splitlines()
    cells=[list(map(int,s.split())) for s in rows[2:42]];j=rows.index('NPOIN= 24')
    xyz=np.array([list(map(float,s.split()[:3])) for s in rows[j+1:j+25]])
    physical={}
    for i,s in enumerate(rows):
        if s.startswith('MARKER_TAG='):
            tag=s.split('=')[1].strip();n=int(rows[i+1].split('=')[1]);physical[tag]=[list(map(int,q.split())) for q in rows[i+2:i+2+n]]
    local={13:[[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]],10:[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]}
    faces={};centers=np.array([xyz[c[1:]].mean(axis=0) for c in cells]);front=centers[:,0]>0
    for i,c in enumerate(cells):
        for idx in local[c[0]]:
            f=[c[j+1] for j in idx];key=tuple(sorted(f));faces.setdefault(key,[]).append((i,f))
    assert all(len(v) in (1,2) for v in faces.values())
    split={k:v for k,v in faces.items() if len(v)==2 and front[v[0][0]]!=front[v[1][0]]}
    assert len(split)>0 and all(np.all(xyz[list(k),0]==0) for k in split)
    details=[]
    for z,selected in enumerate((front,~front)):
        ids=np.flatnonzero(selected);used=np.unique(np.concatenate([cells[i][1:] for i in ids]));mapping={int(n):k for k,n in enumerate(used)}
        tags={name:[f for f in fs if selected[faces[tuple(sorted(f[1:]))][0][0]]] for name,fs in physical.items()}
        interface=[]
        for key,owners in split.items():
            i,f=next(v for v in owners if selected[v[0]]);f=list(f);pts=xyz[f]
            # SU2边界采用指向区域内部的面顺序；与原壳物理壁/远场一致。
            normal=np.cross(pts[1]-pts[0],pts[2]-pts[0])
            if normal@(pts.mean(axis=0)-centers[i])>0:f=f[::-1]
            interface.append([5 if len(f)==3 else 9]+f)
        tags['join%d'%z]=interface
        lines=['NDIME= 3','NELEM= %d'%len(ids)]
        lines += [' '.join(map(str,[cells[i][0]]+[mapping[n] for n in cells[i][1:]])) for i in ids]
        lines+=['NPOIN= %d'%len(used)]+['%.17g %.17g %.17g %d'%(*xyz[n],k) for k,n in enumerate(used)]
        lines+=['NMARK= %d'%len(tags)]
        for name,fs in tags.items():
            lines+=['MARKER_TAG= '+name,'MARKER_ELEMS= %d'%len(fs)]
            lines+=[' '.join(map(str,[f[0]]+[mapping[n] for n in f[1:]])) for f in fs]
        (out/('zone%d.su2'%z)).write_text('\n'.join(lines)+'\n',encoding='utf8')
        details.append(dict(zone=z,motherCells=ids.tolist(),localToMotherNodes=used.tolist(),
          markers={k:len(v) for k,v in tags.items()},meshSHA256=sha(out/('zone%d.su2'%z))))
    assert sum(len(r['motherCells']) for r in details)==40
    save(out/'mesh_mapping.json',dict(originalMeshSHA256=sha(p),regions=details,matchingInterfaceMotherFaces=[list(k) for k in sorted(split)],
      points=24,cells=40,originalWallFaces=8,newGeometryGenerated=False))


def prepare():
    out=BASE/'prepare01';out.mkdir(exist_ok=False);BUILD.mkdir(exist_ok=False)
    single=BASE/'single01';multi=BASE/'multi01';single.mkdir();multi.mkdir()
    split_mesh(multi)
    shutil.copyfile(OLD/'control01/mesh.su2',single/'mesh.su2')
    cfg={k.strip():v.strip() for k,v in [s.split('=',1) for s in (OLD/'control01/case.cfg').read_text(encoding='utf8').splitlines() if '=' in s]}
    # 原数值组合不改；本控制不推进流动。目录与多区元数据只影响真实分派和观测。
    cfg.update(MULTIZONE='NO',MULTIZONE_MESH='NO')
    (single/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    for z in (0,1):
        c=dict(cfg);c.update(MULTIZONE='YES',MESH_FILENAME='zone%d.su2'%z,MARKER_ZONE_INTERFACE='(join0, join1)',
          MARKER_FLUID_INTERFACE='(join%d)'%z,KIND_INTERPOLATION='NEAREST_NEIGHBOR')
        (multi/('zone%d.cfg'%z)).write_text('\n'.join(k+'= '+v for k,v in c.items())+'\n',encoding='utf8')
    driver=dict(SOLVER='MULTIPHYSICS',MATH_PROBLEM='DIRECT',MULTIZONE='YES',MULTIZONE_MESH='NO',
      CONFIG_LIST='(zone0.cfg, zone1.cfg)',MARKER_ZONE_INTERFACE='(join0, join1)',MARKER_FLUID_INTERFACE='(join0, join1)',
      KIND_INTERPOLATION='NEAREST_NEIGHBOR',MULTIZONE_SOLVER='BLOCK_GAUSS_SEIDEL',OUTER_ITER='1',INNER_ITER='1')
    (multi/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in driver.items())+'\n',encoding='utf8')
    cmd=read(OLD/'build/command.json');cmd=[str(CPP) if x.endswith('native\\su2PhaseAwareAir\\compatibility.cpp') else
      x.replace(str(OLD/'build/compatibility.d'),str(BUILD/'control.d')).replace(str(OLD/'build/compatibility.exe'),str(BUILD/'control.exe')) for x in cmd]
    assert str(CPP) in cmd;save(BUILD/'command.json',cmd)
    dep=out/'inputs.d';filtered=[];i=0
    while i<len(cmd):
        x=cmd[i]
        if x in ('-MF','-o'):i+=2;continue
        if x in ('-MD','-Wl,--start-group','-Wl,--end-group') or x.endswith('.a'):i+=1;continue
        filtered.append(x);i+=1
    subprocess.check_call(filtered+['-M','-MT','dep','-MF',str(dep)],cwd=str(OLD_BUILD),timeout=20)
    known=read(GAS/'build02/dependencies.json');oldtree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');bindings={}
    for name in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
        p=Path(name);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
        if p!=CPP:
            if str(p) in known:assert known[str(p)]==h,str(p)
            elif str(p) in oldtree:assert oldtree[str(p)]==h,str(p)
            else:
                rel=p.relative_to(GAS/'source');assert oldtree[str(OLD_SOURCE/rel)]==h,str(p)
        bindings[str(p)]=h
    save(out/'dependencies.json',bindings)
    finish_sources()


def finish_sources():
    """只恢复准确源摘录；已生成配置/网格/依赖不重做、不覆盖失败记录。"""
    out=BASE/'prepare01';known=read(GAS/'build02/dependencies.json')
    oldtree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');bindings=read(out/'dependencies.json')
    for n,h in bindings.items():assert sha(n)==h,n
    source_specs=[(GAS/'source/SU2_CFD/src/variables/CNSVariable.cpp','bool CNSVariable::SetPrimVar',112),
      (OLD_SOURCE/'Common/src/geometry/CGeometry.cpp','void CGeometry::ComputeWallDistance',86),
      (OLD_SOURCE/'SU2_CFD/src/drivers/CMultizoneDriver.cpp','bool CMultizoneDriver::TransferData',50)]
    source_bindings=[];recovered=[]
    for j,(p,needle,count) in enumerate(source_specs):
        h=sha(p);expected=known.get(str(p),oldtree.get(str(p)))
        if expected is None:
            expected=read(GAS/'process_build/report.json')['inputSHA256'][str(p)]
        if h!=expected:
            original=p;archived=ROOT/'outputs/v167_validation/build01/executed_sources'/p.relative_to(OLD_SOURCE)
            assert sha(archived)==expected,str(p)
            recovered.append(dict(currentPath=str(p),currentSHA256=h,executedSourcePath=str(archived),executedSHA256=expected))
            p=archived;h=expected
        rows=p.read_text(encoding='utf8').splitlines();idx=[i for i,s in enumerate(rows) if needle in s];assert len(idx)==1
        a=idx[0];b=min(len(rows),a+count)
        target=out/('%02d_%s.txt'%(j,p.name));excerpt='\n'.join('%6d %s'%(i+1,rows[i]) for i in range(a,b))+'\n'
        if target.exists():assert target.read_text(encoding='utf8')==excerpt
        else:target.write_text(excerpt,encoding='utf8')
        source_bindings.append(dict(path=str(p),sha256=h,firstLine=a+1,lastLine=b))
    save(out/'report.json',dict(passed=True,dependencies=len(bindings),sources=source_bindings,
      recoveredExactSources=recovered,reusedCNSDomainRecovery=True,physicsObjectsRecompiled=False,newDriverSourceSHA256=sha(CPP),flowUpdates=0))


def main(stage):
    if CLOSED:raise SystemExit('v234 control closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda _:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',CPP,DOC/'PLAN.md']
    kind='offline';seconds=120;cwd=BASE
    if stage in ('prepare','prepare2'):
        locks += [OLD/'control01/case.cfg',OLD/'control01/mesh.su2',OLD/'build/command.json',
          GAS/'build02/dependencies.json',ROOT/'outputs/v167_validation/build01/source_tree.json']
        cmd=[sys.executable,str(Path(__file__).resolve()),'task_prepare' if stage=='prepare' else 'task_sources']
    else:
        assert read(BASE/'prepare01/report.json')['passed']
        deps=read(BASE/'prepare01/dependencies.json')
        for n,h in deps.items():assert sha(n)==h,n
        locks += [Path(n) for n in deps]+[BUILD/'command.json']
        for group in read(GAS/'build02/members.json').values():
            for r in group:assert sha(r['path'])==r['sha256'];locks.append(Path(r['path']))
        locks += [GAS/'build02/libSU2core.a',GAS/'build02/libSU2Common.a',GAS/'build02/members.json']
        if stage.startswith('build'):
            kind='build';used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_build*/report.json'));seconds=120-used
            assert seconds>0;cmd=read(BUILD/'command.json');cwd=OLD_BUILD
        elif stage in ('single','multi'):
            assert any(read(p)['passed'] for p in BASE.glob('process_build*/report.json'))
            kind='query';seconds=60-sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='query')
            assert seconds>0;cwd=BASE/(stage+'01');locks += [BUILD/'control.exe']+[p for p in cwd.iterdir() if p.is_file()]
            cmd=[str(BUILD/'control.exe'),'case.cfg','.',stage]
        else:raise ValueError('STAGE_NOT_ADMITTED')
    return guard.supervise(stage,cmd,seconds,128,list(dict.fromkeys(locks)),cwd,kind)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v234 control closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_sources':finish_sources()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
