"""整理自包含研究工程与完整替换文件；隔离单位复现不冒称整箭重算。"""
from pathlib import Path
import json, os, shutil, subprocess, zipfile, hashlib, sys
from screenRocketFullFLeastSquares import ROOT, sha
OUT=ROOT/'outputs/v85_validation';DEL=OUT/'delivery';STAGE=DEL/'isolated/RocketAeroHeat'
CHANGED=['internalFlow/cutCell/buildRocketAugmentedPatchQR.m','internalFlow/cutCell/applyRocketAugmentedPatchQR.m',
 'internalFlow/cutCell/solveRocketAugmentedQRDirection.m','tests/testRocketAugmentedPatchQR.m','tests/studyRocketAugmentedPatchQR.m',
 'runRocketAugmentedQRResearch.m','tests/screenRocketAugmentedPatchQR.py','tests/auditRocketAugmentedPatchQR.py',
 'tests/buildRocketQRWorkbench.py','tests/testRocketQRViewBinding.py','tests/diagnoseRocketQRTemperatureLimit.py',
 'tests/auditRocketWallVelocityWorkbench.py','tests/packageRocketAugmentedQR.py',
 'tests/recordRocketAugmentedQRProgress.py','docs/AUGMENTED_PATCH_QR_V85.md','docs/GITHUB_STAGE_READINESS.md',
 'docs/PROJECT_GOALS.md','docs/validation/v85/progress_augmented_qr.json']


def source_files():
    raw=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT)
    return sorted(ROOT/Path(n.decode('utf8')) for n in raw.split(b'\0') if n and (ROOT/Path(n.decode('utf8'))).is_file())


def copy(file,link=False):
    dest=STAGE/file.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():dest.unlink()
    if link:os.link(file,dest)
    else:shutil.copy2(file,dest)


def stage():
    assert not STAGE.exists();STAGE.mkdir(parents=True)
    base=ROOT/'outputs/v84_validation/delivery/isolated/RocketAeroHeat'
    assert (ROOT/'outputs/v84_validation/delivery/package_verification.json').is_file()
    # 仅复用冻结输出，采用硬链接节约本地空间；源文件另复制，避免改动旧交付。
    for file in (base/'outputs').rglob('*'):
        if file.is_file() and 'isolated_validation' not in file.relative_to(base/'outputs').parts:
            dest=STAGE/file.relative_to(base);dest.parent.mkdir(parents=True,exist_ok=True);os.link(file,dest)
    for file in source_files():copy(file)
    for name in ['outputs/v83_validation/parentlsqrunit02/parent_least_squares_unit_evidence.mat',
                 'outputs/v83_validation/peakweightunit01/peak_weighted_unit_evidence.mat',
                 'outputs/v83_validation/integration01/wall_velocity_integration_evidence.mat',
                 'outputs/v82_validation/exact_krylov01/steady_krylov_evidence.mat',
                 'outputs/v83_validation/krylov26/wall_velocity_krylov_evidence.mat',
                 'outputs/v83_validation/krylov26/wall_velocity_krylov_summary.json']:
        copy(ROOT/name,True)
    for folder in ['qr_screen01','qr_screen02','qr_unit01','qr_real01','rejected01','workbench01','source_epochs']:
        for file in (OUT/folder).rglob('*'):
            if file.is_file():copy(file,True)
    for file in OUT.glob('*.json'):copy(file)
    # 真实当前缓存可由随包Q25/原网格重建，避免重复打包大矩阵。
    for file in OUT.glob('*_stdout.log'):copy(file)
    (STAGE/'AUGMENTED_QR_PACKAGE_README.txt').write_text('''RocketAeroHeat v85 分块QR诊断研究工程

在本目录启动MATLAB，运行：
  main                              原工程入口
  runRocketAugmentedQRResearch       快速单位/稠密参考校核，输出新目录
  runRocketAugmentedQRResearch('outputs/full_qr_reproduction',true)
                                    从随包真实Q25和原网格重建整箭研究，较耗时

当前工作台：outputs/v85_validation/workbench01/index.html。
压力、热流、温度仍是未收敛研究场。线性改善和一步接受未解决非线性温域限制。
正式热流覆盖仍66.5626327424%，新增合格面积0；没有通过新的区域物理或实验验证。
隔离验收范围仅新入口单位/全部局部映射/稠密解复现；没有再隔离重算整箭128步或完整NS接续。
旧历史缓存未全部提供；新整箭入口会重建所需矩阵、因子和当前NS检查。
工作台数据已随包；旧证据里的绝对来源路径用于冻结身份，重新计算请使用上述入口的新目录。
所有原正性、80–2000K、线性和稳态门槛保持。未提交或推送GitHub。
''',encoding='utf8')
    identity={p.relative_to(STAGE).as_posix():sha(p) for p in STAGE.rglob('*') if p.is_file() and p.suffix in ('.m','.js')}
    (DEL/'staged_code_identity.json').write_text(json.dumps(identity,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(stage=str(STAGE),files=sum(p.is_file() for p in STAGE.rglob('*')))))


def verify():
    from scipy.io import loadmat
    import numpy as np
    original=loadmat(OUT/'qr_unit01/augmented_qr_unit_evidence.mat',simplify_cells=True)
    isolated=loadmat(STAGE/'outputs/augmented_qr_reproduction/unit/augmented_qr_unit_evidence.mat',simplify_cells=True)
    comparisons=0
    def equal(a,b):
        nonlocal comparisons
        if isinstance(a,dict):
            assert isinstance(b,dict) and a.keys()==b.keys()
            for k in a:
                if k not in ('elapsedSeconds','constructionSeconds'):equal(a[k],b[k])
        elif isinstance(a,(list,tuple)) or (isinstance(a,np.ndarray) and a.dtype==object):
            assert len(a)==len(b)
            for x,y in zip(a,b):equal(x,y)
        else:
            np.testing.assert_array_equal(a,b);comparisons+=1
    equal(original['cases'],isolated['cases']);equal(original['report'],isolated['report'])
    identity=json.loads((DEL/'staged_code_identity.json').read_text('utf8'))
    assert all(sha(STAGE/name)==h for name,h in identity.items())
    proof=json.loads((STAGE/'outputs/isolated_qr_runtime.json').read_text('utf8'))
    Path(proof['resolvedQR']).resolve().relative_to(STAGE.resolve())
    Path(proof['resolvedEntry']).resolve().relative_to(STAGE.resolve())
    out=dict(status='passed-isolated-QR-entry-and-all-nontiming-small-domain-objects',comparisons=comparisons,
             isolatedProof=proof,allStagedMAndJSSourcesUnchanged=True,fullDomainRecomputedInIsolation=False,
             physicalSolution=False,newProductionSolvedAreaM2=0)
    (DEL/'isolated_comparison.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))


def archive():
    proof=json.loads((DEL/'isolated_comparison.json').read_text('utf8'));assert proof['status'].startswith('passed-isolated-QR')
    identity=json.loads((DEL/'staged_code_identity.json').read_text('utf8'))
    assert all(sha(STAGE/name)==h and (not (ROOT/name).exists() or sha(ROOT/name)==h) for name,h in identity.items() if not name.startswith('outputs/'))
    for file in source_files():
        if file.suffix not in ('.m','.js'):copy(file)
    for file in OUT.glob('*.json'):copy(file)
    archives=[]
    # 隔离验收的临时默认输出不随包，保证解压后缺省入口使用未存在的新目录。
    engineering_files=sorted(p for p in STAGE.rglob('*') if p.is_file() and
        not p.relative_to(STAGE).as_posix().startswith('outputs/augmented_qr_reproduction/'))
    for name,files,base,prefix in [
        ('RocketAeroHeat_v85_research_engineering.zip',engineering_files,STAGE,'RocketAeroHeat/'),
        ('RocketAeroHeat_v85_complete_replacements.zip',[ROOT/n for n in CHANGED],ROOT,'')]:
        dest=DEL/name;assert not dest.exists();entries={}
        with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
            for file in files:
                entry=prefix+file.relative_to(base).as_posix();assert entry not in entries
                z.write(file,entry);entries[entry]=sha(file)
        with zipfile.ZipFile(dest) as z:
            if prefix:
                assert not any(n.startswith(prefix+'outputs/augmented_qr_reproduction/') for n in z.namelist())
                assert prefix+'runRocketAugmentedQRResearch.m' in z.namelist()
            assert z.testzip() is None
            for entry,h in entries.items():
                digest=hashlib.sha256()
                with z.open(entry) as stream:
                    for block in iter(lambda:stream.read(8*1024**2),b''):digest.update(block)
                assert digest.hexdigest()==h,entry
        archives.append(dict(file=str(dest),bytes=dest.stat().st_size,sha256=sha(dest),entries=len(entries),allEntrySHA256=entries))
    out=dict(status='passed-research-engineering-package-isolated-unit-and-all-entry-SHA-CRC',archives=archives,
             defaultRuntimeOutputDirectoryExcluded=True,defaultEntrySourceIdenticalToIsolatedTest=True,
             isolatedReproduction=proof,fullDomainRecomputedInIsolation=False,physicalSolution=False,newProductionSolvedAreaM2=0,githubPushed=False)
    (DEL/'package_verification.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({**out,'archives':[{k:v for k,v in a.items() if k!='allEntrySHA256'} for a in archives]}))


if __name__=='__main__':{'stage':stage,'verify':verify,'archive':archive}[sys.argv[1]]()
