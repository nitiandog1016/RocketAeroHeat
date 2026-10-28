"""自包含交付当前温度投影研究，保留实际证据并隔离复现快速入口。"""
from pathlib import Path
import json,shutil,os,zipfile,hashlib,sys
from screenRocketFullFLeastSquares import ROOT,sha
from packageRocketAugmentedQR import source_files
OUT=ROOT/'outputs/v86_validation';DEL=OUT/'delivery';STAGE=DEL/'isolated/RocketAeroHeat'
CHANGED=['internalFlow/cutCell/cutCellExactEulerBranchJacobian.m','internalFlow/cutCell/buildRocketInterfaceTemperatureRows.m',
 'internalFlow/cutCell/buildRocketTemperatureProjector.m','internalFlow/cutCell/applyRocketTemperatureProjector.m',
 'internalFlow/cutCell/solveRocketTemperatureProjectedDirection.m','tests/testRocketTemperatureProjection.m',
 'tests/studyRocketTemperatureProjectedQR.m','tests/fixtures/rocketExactEulerBranchV85Reference.m',
 'runRocketTemperatureProjectionResearch.m','web/rocket-wall-velocity.js','tests/auditRocketTemperatureProjection.py',
 'tests/plotRocketTemperatureProjection.py','tests/buildRocketTemperatureWorkbench.py','tests/testRocketTemperatureViewBinding.py',
 'tests/diagnoseRocketTemperatureActiveFaces.py','tests/auditRocketWallVelocityWorkbench.py','tests/testRocketWallVelocityWorkbench.js',
 'tests/packageRocketTemperatureProjection.py','tests/recordRocketTemperatureProjectionProgress.py',
 'docs/TEMPERATURE_PROJECTED_QR_V86.md','docs/GITHUB_STAGE_READINESS.md','docs/PROJECT_GOALS.md',
 'docs/validation/v86/progress_temperature_projection.json']


def copy(file,link=False):
    dest=STAGE/file.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():dest.unlink()
    if link:os.link(file,dest)
    else:shutil.copy2(file,dest)


def stage():
    assert not STAGE.exists();STAGE.mkdir(parents=True)
    for file in source_files():copy(file)
    required=['outputs/v83_validation/parentlsqrunit02/parent_least_squares_unit_evidence.mat',
       'outputs/v83_validation/peakweightunit01/peak_weighted_unit_evidence.mat',
       'outputs/v83_validation/riemannreservoirunit02/riemann_reservoir_unit_evidence.mat',
       'outputs/v83_validation/integration01/wall_velocity_integration_evidence.mat',
       'outputs/v82_validation/exact_krylov01/steady_krylov_evidence.mat',
       'outputs/v80_validation/domain02/shell_domain_evidence.mat',
       'outputs/v82_validation/clawpack_exact_euler.py','outputs/v82_validation/probe_exact_riemann.py',
       'outputs/v82_validation/clawpack_Euler.html']
    for case in (26,27,28):
        required.append(f'outputs/v83_validation/krylov{case}/wall_velocity_krylov_evidence.mat')
        required.append(f'outputs/v83_validation/krylov{case}/wall_velocity_krylov_summary.json')
    for name in required:
        file=ROOT/name
        if file.exists():copy(file,True)
        else:assert name.endswith('krylov26/wall_velocity_krylov_summary.json')
    for name in ('temperatureunit86','temperatureunit86b','temperaturefull86','temperatureview86'):
        for suffix in ('_process.json','_stdout.log','.log'):
            copy(ROOT/'outputs/v83_validation'/(name+suffix))
    for name in ('unit02','comparison01','rejected01','workbench01','source_epochs'):
        for file in (OUT/name).rglob('*'):
            if file.is_file():copy(file,True)
    for file in OUT.glob('*'):
        if file.is_file():copy(file)
    (STAGE/'TEMPERATURE_PROJECTION_PACKAGE_README.txt').write_text('''RocketAeroHeat v86 界面温度方向研究工程

在本目录启动MATLAB：
  main
  runRocketTemperatureProjectionResearch
  runRocketTemperatureProjectionResearch('outputs/full_temperature_reproduction',true)

缺省运行真实小域、解析温度行、旧内核多面兼容、稠密受限解校核；输出目录自动使用新编号。
显式true才从随包真实Q26/原网格重建完整当前系统，做自由/受限各128步和原NS有界对照，耗时较长。
旧版本研究缓存并未全部打包；历史入口需要相应历史包。当前新入口的完整输入已提供。
当前工作台：outputs/v86_validation/workbench01/index.html，比较详情在同一工作台内打开。

当前结论：同Q26一次原NS对照，接受步1/16384→1/128，RMS下降0.00609%→0.7848%；不是整箭收敛证明。
受限方向只固定一阶界面温度变化，不添加物理温度边界，不保证非线性温域。
原正性、80–2000K、原NS0.1/内部C1e-3、空间/稳态门槛保持。
当前完整场均为未收敛研究值；正式热流覆盖66.5626327424%，新增合格面积0。
隔离复现仅快速入口与全部非计时小域数值对象，不声称隔离重算整箭128步/NS接续。
科研PNG已查看；无新增浏览器/WebGL视觉验收。完整文件、输入与证据随包，不包含凭据或Git数据。
未提交/推送GitHub；最终校验与运行身份见包旁package_verification.json。
''',encoding='utf8')
    identity={p.relative_to(STAGE).as_posix():sha(p) for p in STAGE.rglob('*') if p.is_file() and p.suffix in ('.m','.js')}
    (DEL/'staged_code_identity.json').write_text(json.dumps(identity,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(stage=str(STAGE),files=sum(p.is_file() for p in STAGE.rglob('*')))))


def verify():
    import numpy as np
    from scipy.io import loadmat
    from scipy.sparse import issparse
    a=loadmat(OUT/'unit02/temperature_projection_unit.mat',simplify_cells=True)
    b=loadmat(STAGE/'outputs/temperature_projection_reproduction/unit/temperature_projection_unit.mat',simplify_cells=True)
    count=0
    def equal(x,y):
        nonlocal count
        if isinstance(x,dict):
            assert isinstance(y,dict) and x.keys()==y.keys()
            for key in x:
                if key not in ('elapsedSeconds','constructionSeconds'):equal(x[key],y[key])
        elif issparse(x):np.testing.assert_array_equal(x.toarray(),y.toarray());count+=1
        elif isinstance(x,(list,tuple)) or (isinstance(x,np.ndarray) and x.dtype==object):
            assert len(x)==len(y)
            for xx,yy in zip(x,y):equal(xx,yy)
        else:np.testing.assert_array_equal(x,y);count+=1
    for key in ('cases','report','Q','gas','scale','rows','fd'):equal(a[key],b[key])
    identity=json.loads((DEL/'staged_code_identity.json').read_text('utf8'));assert all(sha(STAGE/n)==h for n,h in identity.items())
    proof=json.loads((STAGE/'outputs/isolated_temperature_runtime.json').read_text('utf8'))
    for key in ('resolvedEntry','resolvedRows','resolvedBranch'):
        Path(proof[key]).resolve().relative_to(STAGE.resolve())
    out=dict(status='passed-isolated-temperature-entry-and-all-nontiming-small-domain-objects',comparisons=count,
             isolatedProof=proof,allStagedMAndJSSourcesUnchanged=True,fullDomainRecomputedInIsolation=False,
             physicalSolution=False,newProductionSolvedAreaM2=0)
    (DEL/'isolated_comparison.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))


def archive():
    proof=json.loads((DEL/'isolated_comparison.json').read_text('utf8'));assert proof['status'].startswith('passed-isolated-temperature')
    identity=json.loads((DEL/'staged_code_identity.json').read_text('utf8'))
    assert all(sha(STAGE/n)==h and sha(ROOT/n)==h for n,h in identity.items() if not n.startswith('outputs/'))
    for file in source_files():
        if file.suffix not in ('.m','.js'):copy(file)
    for file in OUT.glob('*.json'):copy(file)
    files=sorted(p for p in STAGE.rglob('*') if p.is_file() and not p.relative_to(STAGE).as_posix().startswith('outputs/temperature_projection_reproduction'))
    archives=[]
    for name,entries,base,prefix in [('RocketAeroHeat_v86_research_engineering.zip',files,STAGE,'RocketAeroHeat/'),
                  ('RocketAeroHeat_v86_complete_replacements.zip',[ROOT/n for n in CHANGED],ROOT,'')]:
        dest=DEL/name;assert not dest.exists();identities={}
        with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
            for file in entries:
                member=prefix+file.relative_to(base).as_posix();assert member not in identities;z.write(file,member);identities[member]=sha(file)
        with zipfile.ZipFile(dest) as z:
            assert z.testzip() is None
            if prefix:assert not any(n.startswith(prefix+'outputs/temperature_projection_reproduction') for n in z.namelist())
            for member,h in identities.items():
                digest=hashlib.sha256()
                with z.open(member) as stream:
                    for block in iter(lambda:stream.read(8*1024**2),b''):digest.update(block)
                assert digest.hexdigest()==h
        archives.append(dict(file=str(dest),bytes=dest.stat().st_size,sha256=sha(dest),entries=len(identities),allEntrySHA256=identities))
    out=dict(status='passed-self-contained-temperature-research-package-isolated-unit-and-all-entry-SHA-CRC',archives=archives,
             isolatedReproduction=proof,defaultRuntimeOutputDirectoriesExcluded=True,fullDomainRecomputedInIsolation=False,
             physicalSolution=False,newProductionSolvedAreaM2=0,githubPushed=False)
    (DEL/'package_verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({**out,'archives':[{k:v for k,v in a.items() if k!='allEntrySHA256'} for a in archives]}))


if __name__=='__main__':{'stage':stage,'verify':verify,'archive':archive}[sys.argv[1]]()
