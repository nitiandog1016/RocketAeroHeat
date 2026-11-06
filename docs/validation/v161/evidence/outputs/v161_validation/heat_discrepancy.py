import sys,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from auditSU2IdealAirState import load_state
from su2SSTInterface import wall_geometry
from validateIdealAir import read_rows
from su2ThermalRuntime import sha
base=ROOT/'outputs/v161_validation';case=base/'control'
raw,f,nd=load_state(case/'restart_flow_000000.csv',True,True)
g=wall_geometry(case/'mesh.su2');ids=g['wallIDs'];near=g['normalNeighbor'];d=g['thermalDistance']
T=f['Temperature'];cp=read_rows(case/'ideal_air_audit/oracle.csv')['cp']
lam=f['Laminar_Viscosity'][ids]*cp[ids]/.72
q=lam*(T[near]-300)/d;native=f['Heat_Flux'][ids]
error=abs(q-native)/np.maximum(1,abs(q));ix=np.argsort(error)[-5:][::-1]
rows=[]
for i in ix:
    rows.append(dict(node=int(ids[i]),neighbor=int(near[i]),position=raw[ids[i],1:4].tolist(),neighborTemperature=float(T[near[i]]),temperatureDifference=float(T[near[i]]-300),distance=float(d[i]),conductivity=float(lam[i]),nativeHeat=float(native[i]),recomputedHeat=float(q[i]),relativeError=float(error[i]),oneTemperatureULPHeat=float(np.spacing(T[near[i]])*lam[i]/d[i])))
inputs=[Path(__file__),case/'restart_flow_000000.csv',case/'mesh.su2',case/'ideal_air_audit/oracle.csv']
r=dict(top=rows,failedNodeCount=int(sum(error>=2e-10)),originalThreshold=2e-10,strictPassed=bool(max(error)<2e-10),causeConfirmed=False,inputSHA256={str(p):sha(p) for p in inputs})
(base/'heat_discrepancy.json').write_text(json.dumps(r,indent=2),encoding='utf8');print(json.dumps(rows,indent=2))
