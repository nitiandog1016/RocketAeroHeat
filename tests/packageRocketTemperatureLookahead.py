"""交付v90完整工程与替换文件；复用冻结输入，独立目录只验收小域入口。"""
from pathlib import Path
import json, os, shutil, sys, hashlib, zipfile
import numpy as np
from scipy.io import loadmat
from scipy.sparse import issparse
from packageRocketAugmentedQR import ROOT, source_files
from screenRocketFullFLeastSquares import sha

OUT=ROOT/'outputs/v90_validation'; DEL=OUT/'delivery'; STAGE=DEL/'isolated/RocketAeroHeat'
CHANGED=[
 'runRocketTemperatureLookaheadResearch.m',
 'internalFlow/cutCell/selectRocketTemperatureLookaheadRows.m',
 'internalFlow/cutCell/solveRocketTemperatureActiveDirection.m',
 'internalFlow/cutCell/finishRocketTemperatureActiveDirection.m',
 'tests/testRocketTemperatureLookahead.m','tests/testRocketTemperatureLookaheadIntegration.m',
 'tests/testRocketTemperatureTargetedSelection.m','tests/testRocketTemperatureTargetedSmall.m',
 'tests/studyRocketTemperatureLookahead.m','tests/recoverRocketTemperatureLookaheadStep.m',
 'tests/inspectRocketTemperatureLookaheadRank.m','tests/inspectRocketLookaheadInterfaceRejection.m',
 'tests/fixtures/rocketTemperatureActiveV89Reference.m',
 'tests/auditRocketTemperatureLookahead.py','tests/auditRocketTemperatureTargetedSelection.py',
 'tests/auditRocketTemperatureLookaheadRecovery.py','tests/auditRocketLookaheadInterfaceRejection.py',
 'tests/auditRocketTemperatureProjection.py','tests/diagnoseRocketTemperatureActiveContinuation.py',
 'tests/buildRocketTemperatureActiveWorkbench.py','tests/plotRocketTemperatureActive.py',
 'tests/testRocketWallVelocityWorkbench.js','web/rocket-wall-velocity.js',
 'tests/packageRocketTemperatureLookahead.py','tests/recordRocketTemperatureLookaheadProgress.py',
 'docs/TEMPERATURE_LOOKAHEAD_V90.md','docs/GITHUB_STAGE_READINESS.md','docs/PROJECT_GOALS.md',
 'docs/validation/v90/progress_temperature_lookahead.json']
EXCLUDED=('outputs/temperature_active_reproduction','outputs/temperature_lookahead_reproduction',
 'outputs/temperature_duplicate_reproduction','outputs/temperature_capacity_reproduction','outputs/temperature_breadth_reproduction')


def copy(file,link=False):
    target=STAGE/file.relative_to(ROOT);target.resolve().relative_to(STAGE.resolve());target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists():
        if os.path.samefile(target,file):return
        target.unlink()
    if link:os.link(file,target)
    else:shutil.copy2(file,target)


def stage():
    assert not STAGE.exists();STAGE.mkdir(parents=True)
    old=ROOT/'outputs/v89_validation/delivery';assert json.loads((old/'package_verification.json').read_text('utf8'))['status'].startswith('passed-')
    base=old/'isolated/RocketAeroHeat'
    for file in (base/'outputs').rglob('*'):
        if file.is_file() and not file.relative_to(base).as_posix().startswith(EXCLUDED):
            dest=STAGE/file.relative_to(base);dest.parent.mkdir(parents=True,exist_ok=True);os.link(file,dest)
    for file in source_files():copy(file)
    # 保留真实失败/检查点/恢复/界面诊断；大矩阵仍可由真实初态重建，不重复归档。
    for name in ('unit01','integration03','targeted01','rank01','recovery01','interface01','continuous01',
                 'workbench01','source_epochs','plot_initial01','ui_build_failure01'):
        for file in (OUT/name).rglob('*'):
            if file.is_file() and 'current' not in file.relative_to(OUT/name).parts:copy(file,True)
    manifest=OUT/'continuous01/step01/current/peak_system.json'
    if manifest.exists():copy(manifest,True)
    for file in OUT.iterdir():
        if file.is_file():copy(file)
    for case in ('krylov37','krylov38'):
        for file in (ROOT/'outputs/v83_validation'/case).glob('*'):
            if file.is_file():copy(file,True)
        for file in (ROOT/'outputs/v83_validation').glob(case+'_independent*.json'):copy(file,True)
    for name in ('lookaheadunit90','lookaheadintegration90','lookaheadfull90','lookaheadrank90','targetedunit90','targetedrecovery90','activeview90'):
        for suffix in ('_process.json','_stdout.log','.log'):
            file=ROOT/'outputs/v83_validation'/(name+suffix)
            if file.exists():copy(file,True)
    (STAGE/'TEMPERATURE_LOOKAHEAD_PACKAGE_README.txt').write_text('''RocketAeroHeat v90 前瞻选面与当前状态恢复研究工程

在本目录的MATLAB中运行 main 打开原统一工作流。
runRocketTemperatureLookaheadResearch 默认只复现三组上下温域制造小域，0新完整NS。
runRocketTemperatureActiveResearch 保留旧双温域方向入口。
最新工作台：outputs/v90_validation/workbench01/index.html（离线可打开，亦可用本地静态服务器）。

已实际测试Q37检查点恢复：复用经SHA/坐标/原C-W目标核对的自由128步，新受限128步，
保留124原行及25预测越温单元的所需面，131面中130行激活，原NS缺陷0.0554276通过0.1。
1/8候选仍因内部面97923的2007.55K星区采样拒绝，1/16接受；Q38仍未稳态。
RMS下降6.1871%，正式合格热流覆盖仍66.5626327424%，新增合格面积0。
Tw300K是输入；外缘/恢复/求解壁温未知，不能把研究有限值面积当正式物理解面积。

显式完整研究代码入口：
runRocketTemperatureLookaheadResearch('outputs/my_new_study',true,'outputs/v83_validation/krylov38')
它会重建当前NS/C/W/QR，自由方向和全部当前温度行，不继承旧导数；耗时较长。
本包实际验收默认小域及原真实恢复，未单独重跑上述从Q38重新构造的完整分支。
当前矩阵缓存没有重复归档；原恢复入口和原逐项审计需要重新准备或原工作区缓存。
原失败、来源、诊断和完整中文MATLAB文件均保留。旧历史研究入口可能需要对应旧包。
独立目录验收只证明小域对象与解析路径一致，不冒称整箭隔离重算或新增物理资格。
归档实际身份见包旁package_verification.json。本地小版本未推送GitHub。
''',encoding='utf8')
    identity={p.relative_to(STAGE).as_posix():sha(p) for p in STAGE.rglob('*') if p.is_file() and p.suffix in ('.m','.js')}
    (DEL/'staged_code_identity.json').write_text(json.dumps(identity,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(stage=str(STAGE),files=sum(p.is_file() for p in STAGE.rglob('*')))))


def verify():
    count=0
    def equal(a,b):
        nonlocal count
        if isinstance(a,dict):
            assert isinstance(b,dict) and a.keys()==b.keys()
            for k in a:
                if k not in ('elapsedSeconds','constructionSeconds'):equal(a[k],b[k])
        elif issparse(a):np.testing.assert_array_equal(a.toarray(),b.toarray());count+=1
        elif isinstance(a,(list,tuple)) or isinstance(a,np.ndarray) and a.dtype==object:
            assert len(a)==len(b)
            for x,y in zip(a,b):equal(x,y)
        else:np.testing.assert_array_equal(a,b);count+=1
    a=loadmat(OUT/'targeted01/temperature_targeted_unit.mat',simplify_cells=True)
    b=loadmat(STAGE/'outputs/temperature_lookahead_reproduction/unit/temperature_targeted_small.mat',simplify_cells=True)
    equal(a['cases'],b['cases']);targeted_count=count
    a=loadmat(ROOT/'outputs/v89_validation/active_unit01/temperature_active_unit.mat',simplify_cells=True)
    b=loadmat(STAGE/'outputs/temperature_active_reproduction/unit/temperature_active_unit.mat',simplify_cells=True)
    for k in ('cases','report','bandRows','Q','gas','scale'):equal(a[k],b[k])
    identity=json.loads((DEL/'staged_code_identity.json').read_text('utf8'));assert all(sha(STAGE/k)==h for k,h in identity.items())
    runtime=json.loads((STAGE/'outputs/isolated_lookahead_runtime.json').read_text('utf8'))
    for key,value in runtime.items():
        if key.startswith('resolved'):Path(value).resolve().relative_to(STAGE.resolve())
    out=dict(status='passed-isolated-default-targeted-and-legacy-active-entries-all-nontiming-objects',
        comparisons=count,targetedComparisons=targeted_count,legacyComparisons=count-targeted_count,
        runtime=runtime,allStagedMAndJSSourcesUnchanged=True,fullDomainRecomputedInIsolation=False,
        physicalSolution=False,newProductionSolvedAreaM2=0)
    (DEL/'isolated_comparison.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))


def archive():
    proof=json.loads((DEL/'isolated_comparison.json').read_text('utf8'));assert proof['status'].startswith('passed-isolated-')
    identity=json.loads((DEL/'staged_code_identity.json').read_text('utf8'))
    assert all(sha(STAGE/n)==h and sha(ROOT/n)==h for n,h in identity.items() if not n.startswith('outputs/'))
    for file in source_files():
        if file.suffix not in ('.m','.js'):copy(file)
    for file in OUT.iterdir():
        if file.is_file():copy(file)
    for suffix in ('_process.json','_stdout.log','.log'):
        file=ROOT/'outputs/v83_validation'/('activeisolated90'+suffix)
        if file.exists():copy(file)
    files=sorted(p for p in STAGE.rglob('*') if p.is_file() and not p.relative_to(STAGE).as_posix().startswith(EXCLUDED))
    receipts=[]
    for name,entries,base,prefix in [('research_engineering',files,STAGE,'RocketAeroHeat/'),('complete_replacements',[ROOT/n for n in CHANGED],ROOT,'')]:
        target=DEL/('RocketAeroHeat_v90_'+name+'.zip');assert not target.exists();members={}
        with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
            for file in entries:
                member=prefix+file.relative_to(base).as_posix();assert member not in members
                z.write(file,member);members[member]=sha(file)
        with zipfile.ZipFile(target) as z:
            assert z.testzip() is None
            for member,digest in members.items():
                h=hashlib.sha256()
                with z.open(member) as f:
                    for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
                assert h.hexdigest()==digest
        receipts.append(dict(file=str(target),bytes=target.stat().st_size,sha256=sha(target),entries=len(members),allEntrySHA256=members))
    out=dict(status='passed-v90-engineering-and-complete-replacements-all-entry-CRC-SHA',archives=receipts,
        isolatedReproduction=proof,fullDomainRecomputedInIsolation=False,physicalSolution=False,newProductionSolvedAreaM2=0,githubPushed=False)
    (DEL/'package_verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(status=out['status'],archives=[{k:v for k,v in a.items() if k!='allEntrySHA256'} for a in receipts])))


if __name__=='__main__':{'stage':stage,'verify':verify,'archive':archive}[sys.argv[1]]()
