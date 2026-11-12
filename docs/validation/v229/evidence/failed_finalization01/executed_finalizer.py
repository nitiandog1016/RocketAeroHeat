"""模块职责：归档真实Euler接口证据、保留失败、保护并发，更新当前执行起点。"""
from pathlib import Path
import ast,csv,json,shutil,subprocess,sys,time
import numpy as np
from runPhaseEulerCompatibilityV229 import ROOT,BASE,DOC,OLD_BUILD,read,save,sha,size

OWN=['native/su2PhaseAwareAir/eulerCompatibility.cpp','scripts/runPhaseEulerCompatibilityV229.py',
 'scripts/recoverPhaseEulerObservationV229.py','scripts/su2PhaseEulerContract.py',
 'scripts/verifyPhaseEulerCompatibilityV229.py','scripts/verifyPhaseEulerBuildV229.py',
 'scripts/finalizePhaseEulerNativeV229.py','docs/HYBRID_PHASE_EULER_NATIVE_V229.md',
 'docs/validation/v229/PLAN.md','docs/validation/v229/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def update_docs():
    top='''**2026-11-12 v229当前起点：独立相态Euler已编译，整体原生接口/真实重启通过。** 复用v228四文件/原相态Air，新build02/枚举16/EULER/NONE/5变量；24节点小壳六滑移/四类一般EOS特征、初态/投影pT/非法拒绝通过，独立max7.478e−16，120守恒量第二进程重启逐位同。18身份/3文件拒绝；无流动更新/真实箭体新场/热面积或资格。见docs/HYBRID_PHASE_EULER_NATIVE_V229.md。

首CharacPrimVar仅黏性填写导致Euler读零，首审核差1失败保持；新身份只补已计算边界态观测，不改数值公式，control02通过。恢复构建exit0但输入稳定false：调度main变动；后审核首次期待单文件又失败保持，准确两脚本未调用main/其余AST同及727真实依赖/归档通过，原回执不追改。准备离线23.174秒/构建70.673秒/查询3.250秒，峰1198628864B/1CPU，组清空；OpenMP警告保留。引擎main未单独启动，真实driver及共用归档通过；非完整Jacobian/稳态/热资格。旧源/包/29/engineering54/原黏性资产保持。

**下一次按docs/validation/v229/NEXT_DECISION.md直接唯一真实全局外流有限发展。** 复用build02/准确依赖、v194原64982tet/18027点/19752壁面及v228原面映射/3128面前体责任，另登记固定成熟组合/总步/初0与周期完整保存/时间内存容量/停止合同，当前CLI关闭、无CFD准入。优先原生driver一次初始化/原循环，必要保存/物理逐量残差适配另预算，不重做本轮小壳/物性/鼻部/旧网格/冻结诊断。R_E−300000R_mass先逐节点转换后范数；原−8/稳定/收支/空间外域保持，不续199/制造场、不扫参数或裁态。新原面研究p/T/H/S/U立即接UI，Euler q/固体Tw未知；壁邻—外缘/流管相容后工程热层，强三维收缩底面仍成熟黏性反馈。全面积/科研UI最高、累计后置，小证据GitHub大场E盘，无包/独有删除、不新代理/computer-use，总目标active，下方历史。'''
    p=ROOT/'AGENTS.md';s=p.read_text(encoding='utf8');assert 'v229当前起点' not in s
    first,rest=s.split('\n\n',1);p.write_text(first+'\n\n'+top+'\n\n'+rest,encoding='utf8')
    p=ROOT/'docs/PROJECT_GOALS.md';s=p.read_text(encoding='utf8');parts=s.split('\n\n')
    for i,t in enumerate(parts):
        if t.startswith('**2026-11-12 v228当前起点'):
            parts[i]='**2026-11-12 v229当前起点**：独立相态Euler build02已链接，真实driver的滑移/四类特征、相容能量初态、非法拒绝与24节点120量保存重启通过。首未计算缓存误读、构建调度main输入变动与首次审核失败均保留；准确构建依赖审核另列。没有流动迭代、新热面积或资格。按[v229执行点](validation/v229/NEXT_DECISION.md)直接登记原64982单元唯一有限全局外流；全面积/科研UI并列最高，总目标active。'
        elif t.startswith('当前下一项为**H2独立相态Euler整体原生入口'):
            parts[i]='当前下一项为**H2真实全局外流有限发展及前体连续外缘**。新RAH_PHASE_AIR_EULER已独立构建并通过实际driver整体控制，旧NS/SA身份保持；引擎main尚未单独运行。复用v194原64982tet/18027点/19752壁面与v228映射，先登记固定原生循环、显式0及周期保存/物理能量残差/资源合同。首先服务球鼻之后、收缩之前3128原面/0.317970330711786m²；16.18%仍是目标面积。原黏性32层2微米、全表面责任及收缩/后体/底面的必要黏性反馈保持。'
        elif t.startswith('v224同模型工程热链、v225波系'):
            parts[i]=t.replace('和v228壁列实际回归已实现','、v228壁列实际回归及v229相态Euler整体原生入口已实现').replace('[v228执行点](validation/v228/NEXT_DECISION.md)','[v229执行点](validation/v229/NEXT_DECISION.md)')
        elif t.startswith('每次先读AGENTS.md与本文件第4节'):
            parts[i]='每次先读AGENTS.md与本文件第4节H1—H5，以[v229执行点](validation/v229/NEXT_DECISION.md)覆盖旧下一步：复用build02及原面网格/映射，另登记唯一真实全局外流有限发展，必须保存初0和周期完整场、逐量物理残差及全域/稳定/空间外域合同。完成真实原面研究数据即接UI，审核壁邻—外缘与流管后接3128前体工程层。当前无CFD准入；不重做已过小壳/物性/鼻部/旧网格，不续199/制造场、不裁态或扫参数。全面积/科研UI最高，累计后置。'
    s='\n\n'.join(parts)
    needle='| 原面与体网格 |';i=s.index(needle)
    s=s[:i]+'| 独立相态Euler入口 | v229隔离build02、实际EULER/NONE/滑移/四类特征、120量真实保存重启通过；严格Euler文件合同 | 仅制造控制；原构建输入变化失败保留并单独审核，尚无流动发展/原面外缘或热资格，下一直接原全局外流 |\n'+s[i:]
    p.write_text(s,encoding='utf8')
    p=ROOT/'README.md';s=p.read_text(encoding='utf8');a=s.index('最新[v228');b=s.index('现有工作台“打开结果JSON”继续加载',a)
    new='''最新[v229实际接口验证](docs/HYBRID_PHASE_EULER_NATIVE_V229.md)：相态SU2独立Euler已隔离编译/链接，真实driver的滑移壁、四类特征边界与120量保存重启通过；独立最大缩放差7.478e−16。初次缓存读取和构建输入守卫失败完整保留，恢复结果与适用范围分别记录。**没有流动迭代、新热面积或物理资格。**

下一直接使用原64982单元完整域进行一次有限全局外流计算，再服务3128前体原面/0.317970330711786m²的真实外缘与工程热层。该16.18%是目标面积，已有1306面鼻部热候选保持。按[执行点](docs/validation/v229/NEXT_DECISION.md)先登记实际预算、初0/周期完整保存及物理残差/收支；不重做已过小壳、物性或鼻部。Euler不求壁热流，强三维/收缩/后体/底面仍需成熟黏性反馈。

'''
    p.write_text(s[:a]+new+s[b:],encoding='utf8')

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);evidence=DOC/'evidence';evidence.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    for r in read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(OLD_BUILD/n)==h
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h
    for n,h in read(BASE/'build_audit02/actual_dependencies.json').items():assert sha(n)==h
    assets=read(ROOT/'outputs/v228_validation/final02/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h
    records=[]
    failures={'audit','recover_build','build_audit'}
    for p in BASE.glob('process_*/report.json'):
        stage=p.parent.name[len('process_'):];r=read(p);records.append(r)
        assert not r['liveGroupAfterCleanup'] and r['stopReason'] is None
        assert r['passed']==(stage not in failures),(stage,r['passed'])
        assert r['inputsStable']==(stage!='recover_build')
    assert len(records)==12
    for n in ('audit02','build_audit02','restart_audit01'):assert read(BASE/n/'report.json')['passed']
    # 实际新增观测不改变首控制已有的物理输出：逐位比较，而非仅比较汇总量。
    equivalent={}
    for name in ('initial.csv','wall.csv','saved.csv','restart_flow.csv','refusals.csv'):
        assert (BASE/'control01'/name).read_bytes()==(BASE/'control02'/name).read_bytes(),name
        equivalent[name]=dict(sha256=sha(BASE/'control02'/name),byteIdentical=True)
    def boundary(folder):
        with (BASE/folder/'boundary.csv').open(encoding='utf8') as f:rows=list(csv.DictReader(f))
        return np.array([[float(r[k]) for k in ['node','nx','ny','nz']+['R%d'%j for j in range(5)]] for r in rows)
    assert boundary('control01').tobytes()==boundary('control02').tobytes()
    save(out/'observation_equivalence.json',dict(passed=True,unchangedFields=equivalent,boundaryFluxAndGeometryBitwise=True,
      firstBoundaryPrimitiveCacheUncomputed=True,observedBoundaryStatesOnlyAdded=True))
    update_docs()
    for stage in ('prepare01','recovery01','control01','control02','audit02','restart01','restart_audit01','build_audit02'):
        shutil.copytree(BASE/stage,evidence/stage)
    for stage in ('build01','build02'):
        target=evidence/stage;target.mkdir()
        for p in (BASE/stage).iterdir():
            if p.suffix in ('.json','.d'):shutil.copyfile(p,target/p.name)
    for folder in BASE.glob('process_*'):shutil.copytree(folder,evidence/folder.name)
    shutil.copyfile(out/'observation_equivalence.json',evidence/'observation_equivalence.json')
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        rel=n if n!='.gitattributes' else '.gitattributes.snapshot';q=evidence/'final_source'/rel
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    for n in ('runPhaseEulerCompatibilityV229.py','recoverPhaseEulerObservationV229.py'):
        r=subprocess.run([sys.executable,str(ROOT/'scripts'/n),'control'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
        assert r.returncode!=0 and b'closed' in r.stdout
    totals={kind:sum(r['elapsedSeconds'] for r in records if r['budgetKind']==kind) for kind in ('offline','build','query')}
    assert totals['offline']+time.monotonic()-start<120 and totals['build']<360 and totals['query']<60
    totalbytes=size();assert totalbytes<512*1024**2
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,actualEulerDriverInterfacePassed=True,
      engineMainExecuted=False,flowIterations=0,newHeatAreaM2=0,newQualifiedAreaM2=0,wallHeatFluxSolved=False,
      conservativeRestartValuesBitwise=120,originalBuildGuardFailurePreserved=True,buildAuditAcceptedSeparately=True,
      protectedCount=len(protected),engineeringArtifactsUnchanged=len(engineering['files']),oldObjectsUnchanged=len(old),
      oldDependenciesUnchanged=len(deps),actualNewDependencies=727,assets=assets,
      processSeconds=totals,finalizationSeconds=time.monotonic()-start,outputBytesBeforeReport=totalbytes,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records),allGroupsEmpty=True,
      oldPhysicsAndDefaultUnchanged=True,observationRecoveryOnly=True,MATLABExecuted=False,sharedWebModified=False,
      browserInteractionTested=False,uniqueAssetsDeleted=False,executionEntriesClosed=True)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v229/evidence_manifest.json']+['docs/validation/v229/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),seconds=report['finalizationSeconds'],bytes=totalbytes,processSeconds=totals)))

if __name__=='__main__':main()
