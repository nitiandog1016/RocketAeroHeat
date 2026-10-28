# -*- coding: utf-8 -*-
"""同一真实Newton方向的保守直线与primitive切向路径比较；不将候选当接受解。"""
from pathlib import Path
import hashlib, json, os
import numpy as np
from scipy.io import loadmat
root=Path(__file__).resolve().parents[1]; val=root/'outputs/v82_validation'
source=val/'exact_krylov01/steady_krylov_evidence.mat'; target=val/'primitive_path_comparison.json'
assert not target.exists(); sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest(); before=sha(source)
s=loadmat(source,simplify_cells=True);r=s['result'];gas=s['gas'];cv=gas['R']/(gas['gamma']-1)
trace=r['trace'];trace=[trace] if isinstance(trace,dict) else list(trace);rows=[]
for t in trace:
    q=t['initialQ'];dq=t['linearUpdateU']*r['residualScale'];rho=q[:,0];u=q[:,1:4]/rho[:,None]
    T=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
    dr=dq[:,0];du=(dq[:,1:4]-u*dr[:,None])/rho[:,None]
    dT=(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*dr)/(rho*cv)
    for alpha in [1.,.5]:
        conservative=q+alpha*dq;cr=conservative[:,0];cu=conservative[:,1:4]/cr[:,None]
        ct=(conservative[:,4]-.5*cr*np.sum(cu*cu,axis=1))/(cr*cv)
        pt=T+alpha*dT;pr=rho+alpha*dr
        rows.append(dict(iteration=int(t['iteration']),alpha=alpha,conservativeTemperatureRangeK=[float(ct.min()),float(ct.max())],
                         primitiveTangentTemperatureRangeK=[float(pt.min()),float(pt.max())],primitiveMinDensity=float(pr.min()),
                         primitiveInvalidTemperatureCount=int(np.count_nonzero((pt<gas['temperatureRange'][0])|(pt>gas['temperatureRange'][1]))),
                         primitiveMinimumTemperatureCV=int(pt.argmin())+1,conservativeMinimumTemperatureCV=int(ct.argmin())+1))
assert sha(source)==before
report=dict(status='observed-two-unclipped-paths-of-same-real-Newton-direction-not-full-NS-trials',processID=os.getpid(),rows=rows,
            originalEquationOrBoundaryChanged=False,acceptedSolutionProduced=False,physicalSolution=False,newProductionSolvedAreaM2=0,
            inputSHA256=before,scriptSHA256=sha(__file__))
target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8');print(json.dumps(report))
