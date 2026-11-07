# -*- coding: utf-8 -*-
"""归档本轮原始小证据，绑定最终完整源；不复制大场、不删除独有数据。"""
from pathlib import Path
import hashlib,json,shutil
root=Path(__file__).resolve().parents[2]
src=Path(__file__).resolve().parent
dest=root/'docs/validation/v190/evidence'
assert not dest.exists()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
protected=json.loads((root/'outputs/v176_validation/protected_initial.json').read_text())
assert all(sha(root/p)==h for p,h in protected.items())
parents=[]
for name,code in [('unit01',1),('unit02',1),('unit03',0),('real01',1),('real02',0),('interface01',0),('interface02',0)]:
 p=json.loads((src/(name+'_parent_process.json')).read_text())
 assert p['status']=='terminal' and p['exitCode']==code and p['stoppedByBudget'] is None
 assert p['elapsedSeconds']<p['maximumSeconds']
 parents.append(dict(name=name,**p))
 if not name.startswith('interface'):
  assert json.loads((src/name/'source_after.json').read_text())['allStable']
assert sum(x['elapsedSeconds'] for x in parents)<120
assert json.loads((src/'interface02/source_binding.json').read_text())['allStable']
assert json.loads((src/'real02/report.json').read_text())['passed']
assert json.loads((src/'unit03/report.json').read_text())['passed']
assert sum(p.stat().st_size for p in src.rglob('*') if p.is_file())<16*1024**2
entries=[]
for p in sorted(src.rglob('*')):
 if not p.is_file() or '__pycache__' in p.parts:continue
 rel=p.relative_to(src)
 if p.suffix=='.log':rel=rel.with_suffix('.stdout.txt')
 target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
 assert sha(p)==sha(target)
 entries.append(dict(path=target.relative_to(dest.parent).as_posix(),source=p.relative_to(root).as_posix(),
                     sha256=sha(p),bytes=p.stat().st_size))
source_names=['.gitattributes','AGENTS.md','docs/PROJECT_GOALS.md','airReferenceTemperature.m',
 'calcBoundaryLayerState.m','calcSurfaceHeating.m','evaluateWallConvection.m','calcWallThermalResponse.m',
 'calcFlightHeatingHistory.m','tests/testBoundaryLayerState.m','tests/testAirReferenceTemperature.m',
 'tests/compareReferenceEnthalpyRocket.m','scripts/prepareHybridEulerCase.py','scripts/hybridExternalState.py',
 'tests/testHybridEulerInterface.py','docs/HYBRID_REFERENCE_ENTHALPY_V190.md',
 'docs/validation/v190/PLAN.md','docs/validation/v190/NEXT_DECISION.md']
for name in source_names:
 p=root/name;p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
receipt=dict(schema='RAH.HybridReferenceEnthalpyDelivery.v1',evidence=entries,
 finalSourceSHA256={p:sha(root/p) for p in source_names},
 finalTextNormalizedToLF=True,executedSourceSnapshotsUnchanged=True,
 parents=parents,parentSeconds=sum(p['elapsedSeconds'] for p in parents),budgetSeconds=120,
 independentCSVSeconds=json.loads((src/'independent_csv.json').read_text())['csvVerificationSeconds'],
 protectedFilesUnchanged=len(protected),newPhysicalQualification=False,newSolvedHeatAreaM2=0,
 nativeCFDRuns=0,nativePreparationRuns=0,builds=0,newUIChecks=0,
 adoptedEulerControl='evidence/interface02/control',
 matlabSourceAndInputLocks='five actual runs stable; original raw source snapshots archived',
 pythonSourceBinding='interface01 retrospective copy; interface02 before/after stable',
 simulationScope='synthetic ideal-gas exchange; not production Eilmer reader or model callback',
 runtimeAndLargeFieldsBackedUpOnGitHub=False,uniqueFilesDeleted=0)
(dest.parent/'evidence_manifest.json').write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
print(json.dumps(dict(evidenceFiles=len(entries),evidenceBytes=sum(x['bytes'] for x in entries),
                     parentSeconds=receipt['parentSeconds'],protectedCount=len(protected))))
