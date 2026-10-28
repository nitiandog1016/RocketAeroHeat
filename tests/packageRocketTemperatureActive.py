"""交付当前温域连续研究工程；完整入口输入随包，矩阵缓存可从真实初态重建。"""
from pathlib import Path
import json, os, shutil, sys, zipfile, hashlib, re
from screenRocketFullFLeastSquares import ROOT, sha
from packageRocketAugmentedQR import source_files

OUT = ROOT / 'outputs/v87_validation'
VERSION = 87
DEL = OUT / 'delivery'
STAGE = DEL / 'isolated/RocketAeroHeat'
CHANGED = [
    'internalFlow/cutCell/selectRocketTemperatureBandRows.m',
    'internalFlow/cutCell/solveRocketTemperatureActiveDirection.m',
    'tests/testRocketTemperatureActiveSet.m', 'tests/continueRocketTemperatureActive.m',
    'runRocketTemperatureActiveResearch.m', 'web/rocket-wall-velocity.js',
    'tests/auditRocketTemperatureActive.py', 'tests/buildRocketTemperatureActiveWorkbench.py',
    'tests/plotRocketTemperatureActive.py', 'tests/diagnoseRocketTemperatureActiveContinuation.py',
    'tests/auditRocketWallVelocityWorkbench.py', 'tests/testRocketWallVelocityWorkbench.js',
    'tests/testRocketTemperatureActiveViewBinding.py', 'tests/packageRocketTemperatureActive.py',
    'tests/recordRocketTemperatureActiveProgress.py', 'docs/TEMPERATURE_ACTIVE_CONTINUATION_V87.md',
    'docs/GITHUB_STAGE_READINESS.md', 'docs/PROJECT_GOALS.md',
    'docs/validation/v87/progress_temperature_active.json',
]


def configure(output_folder):
    """选择本轮目录和版本，旧缺省仍指向v87；每包绑定本轮真实首个前驱。"""
    global OUT, DEL, STAGE, CHANGED, VERSION
    folder = Path(output_folder)
    OUT = (folder if folder.is_absolute() else ROOT / folder).resolve()
    OUT.relative_to(ROOT.resolve())
    match = re.fullmatch(r'v(\d+)_validation', OUT.name)
    assert match and OUT.is_dir()
    VERSION = int(match.group(1)); assert VERSION in (87, 88, 89), '请先配置对应版本的目标和交付文件。'
    DEL = OUT / 'delivery'; STAGE = DEL / 'isolated/RocketAeroHeat'
    CHANGED = [n.replace('docs/validation/v87/', f'docs/validation/v{VERSION}/') for n in CHANGED]
    if VERSION == 88:
        CHANGED = [n.replace('docs/TEMPERATURE_ACTIVE_CONTINUATION_V87.md', 'docs/CURRENT_STATE_CONTINUATION_V88.md') for n in CHANGED]
        CHANGED += ['web/viewer.js', 'tests/testViewerFieldState.js', 'tests/buildRocketCoverageIdentityPreview.py',
                    'internalFlow/cutCell/buildRocketTemperatureProjector.m', 'tests/testRocketTemperatureDuplicateRows.m',
                    'tests/fixtures/rocketTemperatureProjectorV87Reference.m',
                    'tests/recoverRocketTemperatureDuplicateStep.m', 'tests/auditRocketTemperatureDuplicateRows.py',
                    'tests/auditRocketTemperatureRecovery.py']
    if VERSION == 89:
        CHANGED = [n.replace('docs/TEMPERATURE_ACTIVE_CONTINUATION_V87.md', 'docs/TEMPERATURE_COMPACT_SUPPORT_V89.md') for n in CHANGED]
        CHANGED += ['internalFlow/cutCell/buildRocketInterfaceTemperatureRows.m',
                    'internalFlow/cutCell/buildRocketTemperatureProjector.m',
                    'internalFlow/cutCell/buildRocketCompactTemperatureProjector.m',
                    'internalFlow/cutCell/applyRocketTemperatureProjector.m',
                    'tests/testRocketTemperatureCapacity.m', 'tests/testRocketTemperatureCompactBreadth.m',
                    'tests/inspectRocketTemperatureCapacity.m', 'tests/auditRocketTemperatureCapacity.py',
                    'tests/inspectRocketTemperatureCompactContext.py',
                    'tests/diagnoseRocketTemperatureTrialCurvature.py',
                    'tests/locateRocketTemperatureStepLimits.py',
                    'tests/rocketTemperatureProjectorReference.py', 'tests/testRocketTemperatureDuplicateRows.m',
                    'tests/fixtures/rocketTemperatureProjectorV87Reference.m']


def copy(file, link=False):
    destination = STAGE / file.relative_to(ROOT)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        assert not os.path.samefile(destination, file), '不能覆盖与源文件相同的硬链接。'
        destination.unlink()
    if link:
        os.link(file, destination)
    else:
        shutil.copy2(file, destination)


def stage():
    assert not STAGE.exists()
    STAGE.mkdir(parents=True)
    for file in source_files():
        copy(file)
    required = [
        'outputs/v83_validation/parentlsqrunit02/parent_least_squares_unit_evidence.mat',
        'outputs/v83_validation/peakweightunit01/peak_weighted_unit_evidence.mat',
        'outputs/v83_validation/riemannreservoirunit02/riemann_reservoir_unit_evidence.mat',
        'outputs/v83_validation/integration01/wall_velocity_integration_evidence.mat',
        'outputs/v82_validation/exact_krylov01/steady_krylov_evidence.mat',
        'outputs/v80_validation/domain02/shell_domain_evidence.mat',
        'outputs/v82_validation/clawpack_exact_euler.py',
        'outputs/v82_validation/probe_exact_riemann.py', 'outputs/v82_validation/clawpack_Euler.html',
    ]
    if VERSION >= 88:
        required.append('outputs/v87_validation/unit01/temperature_active_unit.mat')
    binding = json.loads((OUT / 'view_source_binding.json').read_text('utf8'))
    study = json.loads((OUT / 'continuous01/active_continuation.json').read_text('utf8'))
    initial_case = Path(study['initialStateFile']).parent.name
    cases = list(dict.fromkeys(['krylov28', initial_case] + [x['displayCase'] for x in binding['entries']]))
    for case in cases:
        required += [f'outputs/v83_validation/{case}/wall_velocity_krylov_evidence.mat',
                     f'outputs/v83_validation/{case}/wall_velocity_krylov_summary.json']
        for suffix in ('independent_audit', 'independent_direct_ns_probe', 'independent_riemann_reservoir'):
            file = ROOT / f'outputs/v83_validation/{case}_{suffix}.json'
            if file.exists():
                required.append(str(file.relative_to(ROOT)))
    for name in required:
        copy(ROOT / name, True)
    # 当前入口的完整输入随包。三个大矩阵/解析缓存可由保存的真实端点重新构造，避免重复归档。
    for name in ('unit01', 'active_unit01', 'breadth01', 'capacity01', 'workbench01', 'continuous01', 'source_epochs', 'recovery01', 'duplicate02', 'entry01', 'entry02'):
        for file in (OUT / name).rglob('*'):
            if file.is_file() and not (name == 'continuous01' and 'current' in file.relative_to(OUT / name).parts):
                copy(file, True)
    for file in OUT.iterdir():
        if file.is_file():
            copy(file)
    for name in (f'activeunit{VERSION}', f'activeunit{VERSION}b', f'activefull{VERSION}', f'activerecovery{VERSION}', f'activeview{VERSION}'):
        for suffix in ('_process.json', '_stdout.log', '.log'):
            file = ROOT / 'outputs/v83_validation' / (name + suffix)
            if file.exists():copy(file)
    if VERSION == 89:
        for name in ('capacity89','capacityunit89','capacityfull89'):
            for suffix in ('_process.json','_stdout.log','.log'):
                file=ROOT/'outputs/v83_validation'/(name+suffix)
                if file.exists():copy(file)
        # 保留旧格式诊断的真实失败和修正复测；这些是冻结Q34兼容证据，非v89新流场。
        for name in ('ui_regression_q34','ui_regression_q34b'):
            for file in (OUT/name).rglob('*'):
                if file.is_file() and file.suffix in ('.json','.png') and 'cases' not in file.relative_to(OUT/name).parts:
                    copy(file,True)
    extra=('本轮原第三步曾因精确重复行秩拒绝，失败记录与实际恢复来源分别保留。\n'
           '恢复只复用经当前系统核对的自由128步，新做受限128步；保留全部七行，未截断约束。\n'
           '隔离验收同时运行缺省入口及七行/三条件的重复行测试，不冒称隔离重算整箭。\n'
           'Q34已有17个原10K带内面，超过当前16行选择上限；下一轮先改进表示/内存，不能盲跑。\n'
           if VERSION==88 else
           '本轮保留当前全部温度面，扩展路径按完整QR组支撑存储；448MiB预算和原秩门槛保持。\n'
           '测试另覆盖17/33个独立代数条件，不能冒称真实物理温度边界或新区域资格。\n'
           '隔离默认入口和扩展投影/当前17行复现，不冒称隔离重算整箭三步。\n'
           '附加复现：addpath(\'tests\'); testRocketTemperatureCapacity(\'outputs/my_capacity_check\');\n'
           '          testRocketTemperatureCompactBreadth(\'outputs/my_compact_breadth\');\n'
           if VERSION==89 else '')
    (STAGE / 'TEMPERATURE_ACTIVE_PACKAGE_README.txt').write_text(f'''RocketAeroHeat v{VERSION} 当前温域与真实连续推进研究工程

MATLAB 中进入本工程目录：
  main
  runRocketTemperatureActiveResearch
  runRocketTemperatureActiveResearch('outputs/my_full_active_study',true,'outputs/v83_validation/{initial_case}')

默认只运行双温域小域、自由/活跃方向与稠密最终子空间校核，输出目录自动编号。
显式true及本轮前驱目录才从真实{initial_case}重新建立当前三维系统，最多连续推进三步，耗时较长。
旧两参数调用的缺省前驱Q28保留；对应初态也随包。
当前入口需要的真实初态、原网格和单位输入均随包；旧历史研究入口可能需要对应旧包。
大型当前矩阵缓存未重复归档，可由随包端点重建；已保存独立审计属于原运行证据，
不能在缺失缓存时声称逐项原审计已在隔离目录复演。

工作台：outputs/v{VERSION}_validation/workbench01/index.html。
给定Tw300 K与未知外缘/恢复/求解壁温保持不同身份。云图均为未收敛研究值。
每步重新选择当前温域两端10 K内的全部内部面，仅朝域外方向激活投影。
不添加物理温度边界，不保证非线性候选合法，不宣称不等式最优解。
原NS、正性、80–2000 K、空间/稳态和内部线性门槛保持。
当前正式热流覆盖66.5626327424%，本轮新增合格面积0。
{extra}

已实际检查浏览器三维彩色网格、压力/热流切换、面元探针、算法说明和同窗诊断返回；
对应截图及browser_review.json随包。未查询WebGL错误计数，也未穷举所有浏览器环境。

隔离复现只验收默认新入口和所有非计时小域对象，不冒称整箭三步隔离重算。
完整MATLAB替换文件、源码和对应限制随包；未提交/推送GitHub。
归档身份及验收以包旁package_verification.json为准。
''', encoding='utf8')
    identity = {p.relative_to(STAGE).as_posix(): sha(p) for p in STAGE.rglob('*')
                if p.is_file() and p.suffix in ('.m', '.js')}
    (DEL / 'staged_code_identity.json').write_text(json.dumps(identity, indent=2) + '\n', encoding='utf8')
    print(json.dumps(dict(stage=str(STAGE), files=sum(p.is_file() for p in STAGE.rglob('*')))))


def verify():
    import numpy as np
    from scipy.io import loadmat
    from scipy.sparse import issparse
    original = loadmat(OUT / ('active_unit01/temperature_active_unit.mat' if VERSION==89 else 'unit01/temperature_active_unit.mat'), simplify_cells=True)
    isolated = loadmat(STAGE / 'outputs/temperature_active_reproduction/unit/temperature_active_unit.mat', simplify_cells=True)
    count = 0
    def equal(x, y):
        nonlocal count
        if isinstance(x, dict):
            assert isinstance(y, dict) and x.keys() == y.keys()
            for key in x:
                if key not in ('elapsedSeconds', 'constructionSeconds'):
                    equal(x[key], y[key])
        elif issparse(x):
            np.testing.assert_array_equal(x.toarray(), y.toarray()); count += 1
        elif isinstance(x, (list, tuple)) or (isinstance(x, np.ndarray) and x.dtype == object):
            assert len(x) == len(y)
            for left, right in zip(x, y):
                equal(left, right)
        else:
            np.testing.assert_array_equal(x, y); count += 1
    for key in ('cases', 'report', 'bandRows', 'Q', 'gas', 'scale'):
        equal(original[key], isolated[key])
    active_count=count
    if VERSION == 88:
        original_duplicate=loadmat(OUT/'duplicate02/temperature_duplicate_unit.mat',simplify_cells=True)
        isolated_duplicate=loadmat(STAGE/'outputs/temperature_duplicate_reproduction/temperature_duplicate_unit.mat',simplify_cells=True)
        for key in ('cases','report','duplicates'):equal(original_duplicate[key],isolated_duplicate[key])
    compact_count=0
    if VERSION == 89:
        before_compact=count
        for source, destination in [('unit01','outputs/temperature_capacity_reproduction'),('breadth01','outputs/temperature_breadth_reproduction')]:
            left=loadmat(OUT/source/'temperature_capacity_unit.mat',simplify_cells=True)
            right=loadmat(STAGE/destination/'temperature_capacity_unit.mat',simplify_cells=True)
            for key in ('cases','report'):equal(left[key],right[key])
            if source=='unit01':equal(left['realRows'],right['realRows'])
        compact_count=count-before_compact
    identity = json.loads((DEL / 'staged_code_identity.json').read_text('utf8'))
    assert all(sha(STAGE / n) == h for n, h in identity.items())
    runtime = json.loads((STAGE / 'outputs/isolated_active_runtime.json').read_text('utf8'))
    keys=('resolvedEntry','resolvedActive','resolvedRows')+ (('resolvedDuplicateTest','resolvedProjector') if VERSION==88 else ())+ (('resolvedCapacity','resolvedCompact','resolvedBreadth') if VERSION==89 else ())
    for key in keys:
        Path(runtime[key]).resolve().relative_to(STAGE.resolve())
    proof = dict(status='passed-isolated-active-temperature-entry-and-all-nontiming-small-domain-objects',
                 comparisons=count,activeUnitComparisons=active_count,duplicateUnitComparisons=count-active_count-compact_count,
                 compactUnitComparisons=compact_count,
                 isolatedProof=runtime, allStagedMAndJSSourcesUnchanged=True,
                 fullDomainRecomputedInIsolation=False, physicalSolution=False, newProductionSolvedAreaM2=0)
    (DEL / 'isolated_comparison.json').write_text(json.dumps(proof, indent=2) + '\n', encoding='utf8')
    print(json.dumps(proof))


def archive():
    proof = json.loads((DEL / 'isolated_comparison.json').read_text('utf8'))
    assert proof['status'].startswith('passed-isolated-active')
    identity = json.loads((DEL / 'staged_code_identity.json').read_text('utf8'))
    assert all(sha(STAGE / n) == h and sha(ROOT / n) == h for n, h in identity.items() if not n.startswith('outputs/'))
    for file in source_files():
        if file.suffix not in ('.m', '.js'):
            copy(file)
    for file in OUT.glob('*.json'):
        copy(file)
    excluded=('outputs/temperature_active_reproduction','outputs/temperature_duplicate_reproduction',
              'outputs/temperature_capacity_reproduction','outputs/temperature_breadth_reproduction')
    files = sorted(p for p in STAGE.rglob('*') if p.is_file()
                   and not p.relative_to(STAGE).as_posix().startswith(excluded))
    receipts = []
    for name, entries, base, prefix in [
        (f'RocketAeroHeat_v{VERSION}_research_engineering.zip', files, STAGE, 'RocketAeroHeat/'),
        (f'RocketAeroHeat_v{VERSION}_complete_replacements.zip', [ROOT / n for n in CHANGED], ROOT, ''),
    ]:
        target = DEL / name; assert not target.exists(); members = {}
        with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED, compresslevel=3, allowZip64=True) as archive:
            for file in entries:
                member = prefix + file.relative_to(base).as_posix(); assert member not in members
                archive.write(file, member); members[member] = sha(file)
        with zipfile.ZipFile(target) as archive:
            assert archive.testzip() is None
            assert not any(n.startswith(tuple(prefix+x for x in excluded)) for n in archive.namelist())
            for member, digest in members.items():
                check = hashlib.sha256()
                with archive.open(member) as stream:
                    for block in iter(lambda: stream.read(8 * 1024**2), b''):
                        check.update(block)
                assert check.hexdigest() == digest
        receipts.append(dict(file=str(target), bytes=target.stat().st_size, sha256=sha(target),
                             entries=len(members), allEntrySHA256=members))
    out = dict(status='passed-self-contained-current-temperature-package-isolated-unit-and-all-entry-SHA-CRC',
               archives=receipts, isolatedReproduction=proof, defaultRuntimeOutputDirectoriesExcluded=True,
               fullDomainRecomputedInIsolation=False, physicalSolution=False, newProductionSolvedAreaM2=0, githubPushed=False)
    (DEL / 'package_verification.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
    print(json.dumps({**out, 'archives': [{k: v for k, v in a.items() if k != 'allEntrySHA256'} for a in receipts]}))


if __name__ == '__main__':
    if len(sys.argv) > 2:
        configure(sys.argv[2])
    {'stage': stage, 'verify': verify, 'archive': archive}[sys.argv[1]]()
