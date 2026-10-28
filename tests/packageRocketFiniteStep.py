"""生成自包含研究工程、隔离运行目录及完整替换文件；归档不等于物理验证。"""
from pathlib import Path
import json,shutil,subprocess,zipfile,hashlib,sys
from recordRocketFiniteStepProgress import SOURCES,sha
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v84_validation';DELIVERY=OUT/'delivery';STAGE=DELIVERY/'isolated'/'RocketAeroHeat'

def sources():
 names=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT).split(b'\0')
 files={ROOT/Path(x.decode('utf8')) for x in names if x};return sorted(p for p in files if p.is_file())
def copy(p):
 target=STAGE/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
def stage():
 assert not STAGE.exists();STAGE.mkdir(parents=True)
 for p in sources():copy(p)
 required=['outputs/v83_validation/krylov24/wall_velocity_krylov_evidence.mat','outputs/v83_validation/krylov25/wall_velocity_krylov_evidence.mat','outputs/v80_validation/domain02/shell_domain_evidence.mat',
  'outputs/v82_validation/clawpack_exact_euler.py','outputs/v82_validation/probe_exact_riemann.py','outputs/v82_validation/clawpack_Euler.html']
 for name in required:copy(ROOT/name)
 for name in ['diagnosis01','patch01','patch02','patch03','conditioning01','unit01','returned_model01','workbench01']:
  for p in (OUT/name).rglob('*'):
   if p.is_file():copy(p)
 for p in OUT.glob('*.json'):copy(p)
 for ext in ['png','pdf']:copy(OUT/f'finite_step_comparison.{ext}')
 for case in ['full01','full02']:
  for p in (OUT/case).rglob('*'):
   if not p.is_file():continue
   # 保留真实接受状态、探针和选择证据；可重建的大矩阵/局部因子不重复打包。
   relative=p.relative_to(OUT/case)
   keep=p.suffix=='.json' or any(part.startswith(('krylov','step')) for part in relative.parts) or p.name=='full_comparison.mat'
   if keep:copy(p)
 for p in OUT.glob('source_epochs/*'):
  if p.is_file():copy(p)
 readme='''RocketAeroHeat v84 有限步诊断研究工程（本地交付）

在解压后的 RocketAeroHeat 目录启动 MATLAB：
  main                            原工程入口
  runRocketFiniteStepResearch     自包含真实几何六步小复现及代数检查
  runRocketFullMinimumStateComparison  从随包 Q24 和原网格重建完整整箭对照，运行较久

outputs/v84_validation/workbench01/index.html 为五步骤科研工作台。
完整源代码、测试、模型、原 Q24 检查点和网格均随包提供。无需旧机缓存矩阵；整箭入口会重新构建。
所有可运行入口使用当前工程目录，计算输出须写入新目录，勿覆盖旧证据或重复启动正在运行的任务。

研究状态：局部六步改善明显，但后续仍停滞；整箭两组只接受一步、第二方向即拒绝，未满足连续多步改善。
研究方向选择不替换默认求解器。详细终态见 docs/FINITE_STEP_DIAGNOSIS_V84.md。
内部C目标1e-3、原NS方向0.1、正性、80–2000K和原稳态/守恒门槛保持。
未收敛或未通过物理验证的场量仍是研究值，正式面积新增0。不能将包可运行视为物理解验证。
隔离运行仅复核小复现和完整Q24基态；完整多步对照使用本项目本轮实际记录，未另作整轮隔离重复计算。
历史所有版本的缓存并未随包提供；本轮新入口所需输入已包含。
未提交或推送GitHub。ZIP全文件校验见包旁 package_verification.json。
'''
 (STAGE/'FINITE_STEP_PACKAGE_README.txt').write_text(readme,encoding='utf8')
 tested_sources={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sources() if p.suffix in ['.m','.js']}
 (DELIVERY/'staged_code_identity.json').write_text(json.dumps(tested_sources,indent=2),encoding='utf8')
 print(json.dumps({'stage':str(STAGE),'files':sum(p.is_file() for p in STAGE.rglob('*'))}))

def archive():
 proof=STAGE/'outputs/isolated_validation/isolated_package.json';assert proof.exists()
 p=json.loads(proof.read_text('utf8'));assert p['status'].startswith('passed-isolated-package')
 independent=json.loads((DELIVERY/'isolated_comparison.json').read_text('utf8'))
 assert independent['status'].startswith('passed-isolated-six-step') and independent['isolatedProofSHA256']==sha(proof)
 tested=json.loads((DELIVERY/'staged_code_identity.json').read_text())
 for name,digest in tested.items():assert sha(ROOT/name)==digest and sha(STAGE/name)==digest, name
 # 仅同步运行后形成的最终记录/文档；已测试MATLAB和前端源不能悄悄更换。
 for file in sources():
  if file.suffix not in ['.m','.js']:copy(file)
 for file in OUT.glob('*.json'):copy(file)
 shutil.copy2(proof,DELIVERY/'isolated_package.json')
 archives=[]
 for filename,files,prefix in [
  ('RocketAeroHeat_v84_research_engineering.zip',sorted(p for p in STAGE.rglob('*') if p.is_file()),'RocketAeroHeat/'),
  ('RocketAeroHeat_v84_complete_replacements.zip',[ROOT/name for name in SOURCES if (ROOT/name).exists()]+[ROOT/'docs/validation/v84/progress_finite_step.json'],'')]:
  dest=DELIVERY/filename;assert not dest.exists();entries={}
  with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
   for file in files:
    name=prefix+file.relative_to(STAGE if prefix else ROOT).as_posix();assert name not in entries
    z.write(file,name);entries[name]=sha(file)
  with zipfile.ZipFile(dest) as z:
   assert z.testzip() is None
   for name,expected in entries.items():
    h=hashlib.sha256()
    with z.open(name) as f:
     for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    assert h.hexdigest()==expected,name
  archives.append({'file':str(dest.relative_to(ROOT)),'bytes':dest.stat().st_size,'sha256':sha(dest),'entries':len(entries),'allEntrySHA256':entries})
 report={'status':'passed-self-contained-engineering-package-isolated-local-run-and-all-entry-verification',
  'archives':archives,'isolatedMATLAB':p,'isolatedIndependentComparison':independent,'fullMultistepRecomputedInIsolation':False,'newProductionSolvedAreaM2':0,'githubPushed':False}
 (DELIVERY/'package_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
 print(json.dumps({'status':report['status'],'archives':[{k:v for k,v in a.items() if k!='allEntrySHA256'} for a in archives]}))
if __name__=='__main__':stage() if sys.argv[1]=='stage' else archive()
