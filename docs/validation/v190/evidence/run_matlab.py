from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'scripts'))
from validateIdealAir import run
folder=Path(__file__).resolve().parent
name=sys.argv[1]
assert name in ('unit01','unit02','unit03','real01','real02')
case=folder/name;case.mkdir(exist_ok=False)
files=list(root.glob('*.m'))+list((root/'tests').glob('*.m'))
before={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(case/'source_before.json').write_text(json.dumps(before,indent=2)+'\n')
for p in [root/'airReferenceTemperature.m',root/'calcBoundaryLayerState.m',root/'calcSurfaceHeating.m',root/'evaluateWallConvection.m',root/'calcWallThermalResponse.m',root/'calcFlightHeatingHistory.m',root/'tests/testBoundaryLayerState.m',root/'tests/testAirReferenceTemperature.m']:
 d=case/'executed_source'/p.relative_to(root);d.parent.mkdir(parents=True,exist_ok=True);d.write_bytes(p.read_bytes())
expression="restoredefaultpath;maxNumCompThreads(1);addpath(pwd);addpath(fullfile(pwd,'tests'));r=testAirReferenceTemperature();fid=fopen('outputs/v190_validation/"+name+"/report.json','w');fprintf(fid,'%s',jsonencode(r));fclose(fid);"
if name.startswith('real'):
 p=root/'outputs/v104_validation/baseline02/current_case.mat';before[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
 (case/'source_before.json').write_text(json.dumps(before,indent=2)+'\n')
 p=root/'tests/compareReferenceEnthalpyRocket.m';(case/'executed_source/tests'/p.name).write_bytes(p.read_bytes())
 expression=expression.replace('r=testAirReferenceTemperature();',"r=compareReferenceEnthalpyRocket('outputs/v190_validation/"+name+"');")
try:
 run(folder,name+'_parent',['D:/matlab/bin/matlab.exe','-batch',expression],40 if name=='real02' else 60)
finally:
 stable={p:hashlib.sha256((root/p).read_bytes()).hexdigest()==h for p,h in before.items()}
 (case/'source_after.json').write_text(json.dumps(dict(allStable=all(stable.values()),changed=[p for p,ok in stable.items() if not ok]),indent=2)+'\n')
