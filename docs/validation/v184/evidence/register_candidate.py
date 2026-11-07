"""事前登记唯一真实域显式候选，不准备或启动求解器。"""
from pathlib import Path
import copy,json,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from eilmerExplicitAdmission import sha,physical_hash,audit_admission
from prepareEilmerRocketCase import validate_condition
base=Path(__file__).resolve().parent;candidate=base/'candidate.json';assert not candidate.exists()
source=ROOT/'outputs/v178_validation/real01';m=json.loads((source/'case_manifest.json').read_text())
c=copy.deepcopy(m['condition']);c.pop('maxNewtonSteps')
evidence={
    'sourceManifest':source/'case_manifest.json',
    'sourcePreparation':ROOT/'outputs/v178_validation/real01_prepare/report.json',
    'sourceGeometry':ROOT/'outputs/v178_validation/real01_geometry.json',
    'explicitInterface':ROOT/'outputs/v183_validation/verify01/report.json',
    'explicitRun':ROOT/'outputs/v183_validation/small01_recovery/report.json'}
r=dict(schema='RAH.EilmerExplicitCandidate.v1',registeredUTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
       profile='rah-explicit-local-development-v1',case=str(base/'real01'),sourceCase=str(source),runRecord=str(base/'real01_run'),
       maxSteps=200,checkpointStep=100,sourceMeshSHA256=m['mesh']['sourceSHA256'],physicalConditionSHA256=physical_hash(c),
       evidence={k:dict(path=str(p),sha256=sha(p)) for k,p in evidence.items()},
       budget=dict(nativeWallClockSeconds=600,runEntrySeconds=660,parentSeconds=720,prepareSeconds=120,
                   cpuCount=2,memoryMiB=8192,outputMiB=1024,automaticContinuation=0),
       screen=dict(window=20,minimumIndependentSamples=40,minMassAndEnergyMeanLogImprovement=0.05),
       physicalQualification=False,newQualifiedAreaM2=0)
candidate.write_text(json.dumps(r,indent=2)+'\n')
c['explicit']=dict(profile=r['profile'],maxSteps=200,checkpointStep=100,candidate=str(candidate),candidateSHA256=sha(candidate))
c=validate_condition(c);audit_admission(dict(condition=c,mesh=m['mesh']),base/'real01')
(base/'condition.json').write_text(json.dumps(c,indent=2)+'\n')
print(json.dumps(dict(candidate=str(candidate),sha256=sha(candidate),nativeRuns=0)))
