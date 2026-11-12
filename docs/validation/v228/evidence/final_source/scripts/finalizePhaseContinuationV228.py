"""模块职责：归档实际壁列回归和独立Euler源码准备，更新当前目标并保护并发。"""
from pathlib import Path
import json,ast,shutil,subprocess,sys,time,os
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from runPhaseContinuationV228 import SOURCES
BASE=ROOT/'outputs/v228_validation';DOC=ROOT/'docs/validation/v228'
OWN=['tests/testPhaseNoseRecoveryV228.m','scripts/runPhaseContinuationV228.py','scripts/preparePhaseEulerV228.py',
 'scripts/stagePhaseEulerSourceV228.py','scripts/finalizePhaseContinuationV228.py',
 'docs/HYBRID_PHASE_CONTINUATION_V228.md','docs/validation/v228/PLAN.md','docs/validation/v228/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def update_docs():
    # 用最新执行点替换旧进度文字，保留H1—H5和全部长期目标，不继续堆叠过期“下一步”。
    p=ROOT/'docs/PROJECT_GOALS.md';s=p.read_text(encoding='utf8');parts=s.split('\n\n')
    for i,text in enumerate(parts):
        if text.startswith('**2026-11-12 v227最新覆盖说明'):
            parts[i]='**2026-11-12 v228当前起点**：v227壁列修复已通过唯一30°BVP实际回归，主工作台15字段/两区域/全19752面数值一致；原失败保留。已有鼻部1306面/0.1790893%热候选保持，本轮新面积/物理资格0。下一直接实施同成熟SU2的独立相态Euler入口，四文件隔离适配、原完整域提案和精确映射已准备，尚未编译/运行。按[v228执行点](validation/v228/NEXT_DECISION.md)先登记隔离构建与唯一整体原生兼容预算，再进入真实全局外流与3128面前体相容入口；不重做鼻部/物性/旧网格或盲续199。总目标active。'
        elif text.startswith('2026-11-07，用户明确'):
            parts[i]='2026-11-07，用户明确将主技术路线调整为**低成本全局外流＋工程边界层＋局部黏性CFD**。**单工况全面积气动热与统一科研UI并列最高。** 自研集中于相容接口/分区/映射/独立验证，CFD复用成熟内核。累计热量、轨迹和跨流态后置，原几何/SI/各类温度身份、近壁资产与所有科学门槛保持。'
        elif text.startswith('当前下一项是 **H2'):
            parts[i]='当前下一项为**H2独立相态Euler整体原生入口及前体连续外缘**。旧v220 RAH_PHASE_AIR仅准NS/SA，新RAH_PHASE_AIR_EULER必须隔离注册、构建并整体实测；不能只改SOLVER或借用旧NS资格。复用v194原64982tet/18027点/19752壁面的无黏网格，不改黏性32层2微米资产。先接入拟合球鼻之后、收缩之前3128原面/0.317970330711786m²；该16.18%是目标面积，不是新增解。其余原面全部保留责任，后体/底面/强三维继续成熟黏性反馈。'
        elif text.startswith('v224同模型热链、v225波系'):
            parts[i]='v224同模型工程热链、v225波系、v226鼻部外缘、v227真实局部相似热剖面和v228壁列实际回归已实现。下一直接通过独立相态Euler提供全局压力/熵/总焓/三维方向，审核壁邻—外缘、相邻标签尺度及完整入口后接前体工程热层。局部相似剖面不自动取得CFD入口资格，不能整后体固定驻点熵或以独立切楔拼全三维。详见[v228执行点](validation/v228/NEXT_DECISION.md)。'
        elif text.startswith('每次先读AGENTS.md与本文件第4节'):
            parts[i]='每次先读AGENTS.md与本文件第4节H1—H5，以[v228执行点](validation/v228/NEXT_DECISION.md)覆盖旧下一步：复用四文件隔离适配及准确源，另预算一次构建/整体Euler滑移、原始量、一般EOS混合特征、同能量初场保存重启与非法拒绝控制。通过才原64982全域有限发展→3128前体外缘/流管/工程边界层。当前无构建/原生/CFD准入，不重做30°/31站/物性/波系/旧网格/冻结诊断；不续199、不裁状态或品牌轮换。全面积/科研UI最高，累计后置。'
    s='\n\n'.join(parts).replace('当前混合链路尚未实现或获物理资格','当前混合链路仅部分模块已实现，尚未完整闭合或获物理资格');p.write_text(s,encoding='utf8')
    p=ROOT/'README.md';s=p.read_text(encoding='utf8')
    if '最新[v228执行与准备]' in s:return
    start=s.index('最新为[v227');end=s.index('现有工作台通过“打开结果JSON”读取[v223',start)
    replacement='''最新[v228执行与准备](docs/HYBRID_PHASE_CONTINUATION_V228.md)：v227壁列修复的实际30°BVP和主工作台15字段/区域面积通过，原失败保留。已有1306面鼻部热候选/0.1790893%/866.146994W保持，本轮新面积及物理资格0。

下一优先接入3128个前体原面/0.317970330711786m²，约占全表面16.18%；这是目标面积。相态SU2独立Euler四文件补丁、64982tet原完整域配置提案与19752面精确映射已准备，**尚未编译/运行**，旧NS引擎不可直接启动。按[执行点](docs/validation/v228/NEXT_DECISION.md)直接隔离构建/整体兼容控制，不重做物性、鼻部或旧网格；必要黏性反馈仍负责收缩/后体/底面。

现有工作台“打开结果JSON”继续加载[鼻部15字段](docs/validation/v227/evidence/recovered/rocket_workbench.json)：本轮实际主导出的全部数值/状态与它一致。q、p、温度、剪应力、Cf、真实剖面积分及来源分开显示；给定Tw和未知固体Tw不混合。共享WebGL未修改，无本轮浏览器交互验收。

小版本同步公开GitHub；完整源码直接更新，不默认生成ZIP/工程包。大体场及重复JSON留E盘，必要紧凑证据公开，保护并发和独有数据。单时刻全面积与科研UI最高，累计热量后置，总目标active。

'''
    p.write_text(s[:start]+replacement+s[end:],encoding='utf8')

def main():
    start=time.monotonic();out=BASE/'final02';out.mkdir(exist_ok=False);evidence=DOC/'evidence';evidence.mkdir(exist_ok=True)
    assert not list(evidence.iterdir())
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json');assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    candidate=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidate['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    assets=read(ROOT/'outputs/v227_validation/final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    for n in SOURCES:assert sha(ROOT/n)==sha(BASE/'matlab01/executed_source'/n),n
    records=[]
    for p in BASE.glob('process_*/report.json'):
        r=read(p);assert r['inputsStable'] and not r['liveGroupAfterCleanup'];records.append(r)
    for stage,passed in [('matlab01',True),('prepare01',False),('prepare02',False),('prepare03',True),('audit01',True),('stage01',True)]:
        assert read(BASE/('process_'+stage)/'report.json')['passed']==passed
    for stage in ('matlab01','prepare03','audit01','stage01'):assert read(BASE/stage/'report.json')['passed']
    update_docs()
    # 大型重复UI JSON保留本地，逐值对照已有公开版本后只公开等价证据。
    actual=read(BASE/'matlab01/workbench.json');prior=read(ROOT/'docs/validation/v227/evidence/recovered/rocket_workbench.json')
    for a,b in zip(actual['fields'],prior['fields']):
        assert a['key']==b['key'] and a['values']==b['values'] and a['status']==b['status']
    assert actual['vertices']==prior['vertices'] and actual['faces']==prior['faces'] and actual['region']==prior['region']
    # 旧纯函数测试采用固定文件名；同盘硬链接只提供名字，不复制或重算字段。
    os.link(str(BASE/'matlab01/workbench.json'),str(out/'rocket_workbench.json'))
    assert os.path.samefile(str(BASE/'matlab01/workbench.json'),str(out/'rocket_workbench.json'))
    run=subprocess.run(['node',str(ROOT/'tests/testPhaseNoseHeatingViewerV227.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend.txt').write_bytes(run.stdout);assert run.returncode==0,run.stdout
    save(out/'viewer_equivalence.json',dict(passed=True,fields=15,originalFaces=19752,numericValuesAndStatusIdentical=True,
      actualExporterExecuted=True,sourceSHA256={n:sha(ROOT/n) for n in SOURCES},browserInteractionTested=False,
      sharedViewerSHA256=sha(ROOT/'web/viewer.js'),testSHA256=sha(ROOT/'tests/testPhaseNoseHeatingViewerV227.js'),
      localActualJSON=dict(path=str(BASE/'matlab01/workbench.json'),sha256=sha(BASE/'matlab01/workbench.json'),bytes=(BASE/'matlab01/workbench.json').stat().st_size),
      publicEquivalentJSON='docs/validation/v227/evidence/recovered/rocket_workbench.json'))
    for folder in ('matlab01','prepare01','prepare02','prepare03','audit01','stage01'):
        target=evidence/folder;target.mkdir()
        for p in (BASE/folder).iterdir():
            if folder=='matlab01' and p.name=='workbench.json':continue
            if folder=='stage01' and p.name=='isolated_source':continue
            if p.is_dir():shutil.copytree(p,target/p.name)
            else:shutil.copyfile(p,target/p.name)
    for p in BASE.glob('process_*'):
        if p.is_dir():shutil.copytree(p,evidence/p.name)
    raw=BASE/'process_matlab01/stdout.txt';decoded=raw.read_bytes().decode('gb18030',errors='replace')
    (evidence/'process_matlab01/stdout_utf8.txt').write_text(decoded,encoding='utf8')
    save(evidence/'process_matlab01/decoding.json',dict(encodingAssumption='gb18030',replacementCharacters=decoded.count('\ufffd'),originalBytesRetained=True))
    for p in ('frontend.txt','viewer_equivalence.json'):shutil.copyfile(out/p,evidence/p)
    shutil.copytree(BASE/'final01',evidence/'failed_finalization01')
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        rel=Path(n) if n!='.gitattributes' else Path('.gitattributes.snapshot')
        q=evidence/'final_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    # 保留新补丁实际生成源与准确依赖的字节来源，不为归档再执行数值或引擎。
    patchManifest=read(BASE/'stage01/manifest.json')
    for rel,v in patchManifest.items():assert sha(v['isolatedPath'])==v['patchedSHA256'] and sha(v['original']['path'])==v['original']['sha256']
    closed=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseContinuationV228.py'),'matlab01'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert closed.returncode!=0 and b'closed' in closed.stdout
    total={k:sum(r['elapsedSeconds'] for r in records if r['budgetKind']==k) for k in ('matlab','offline')}
    assert total['matlab']<150 and total['offline']<120
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<128*1024**2
    elapsed=time.monotonic()-start;assert elapsed<120
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,newHeatAreaM2=0,newQualifiedAreaM2=0,
      currentMatlabWallCorrectionExecuted=True,mainExporterExecuted=True,all31StationsRerun=False,
      protectedCount=29,sharedEngineeringArtifactsUnchanged=len(engineering['files']),oldObjectsStable=len(old),dependenciesStable=len(deps),
      assets=assets,firstFailuresRetained=True,sourceAdapterStagedNotCompiled=True,oldNSGuardUnchanged=True,
      target=read(BASE/'prepare03/report.json'),processBudgetSeconds=total,finalizationSeconds=elapsed,
      outputBytesBeforeFinalReport=size,outputByteCountIncludesHardlinkAlias=True,CFDRuns=0,nativeQueries=0,builds=0,validationEntryClosed=True,
      sharedWebSourceCopied=False,browserInteractionTested=False,uniqueAssetsDeleted=False)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v228/evidence_manifest.json']+['docs/validation/v228/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,files=len(manifest),seconds=elapsed,outputBytes=size,processBudgetSeconds=total)))
if __name__=='__main__':main()
