# -*- coding: utf-8 -*-
"""独立复算已拒绝界面的理想气体Riemann压力根与波系；不接受越域候选。"""
from pathlib import Path
import ast, hashlib, json, math, os, re, sys
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq, bisect

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';REF=ROOT/'outputs/v82_validation'
case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_interface.json');assert not target.exists()
file=VAL/case/'rejected_interface_evidence.mat';curves=REF/'clawpack_exact_euler.py';probe=REF/'probe_exact_riemann.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
bindings={str(p.relative_to(ROOT)):sha(p) for p in [file,curves,probe]}
assert sha(curves)=='2a926d59d4935807b93c425aa617bf64d64858e65231fd556b6ae6202d0ec2ef'
allowed={'pospart','sound_speed','beta','integral_curve_1','integral_curve_3','hugoniot_locus_1','hugoniot_locus_3'}
nodes=[n for n in ast.parse(curves.read_text(encoding='utf8')).body if isinstance(n,ast.FunctionDef) and n.name in allowed]
assert {n.name for n in nodes}==allowed
for definition in nodes:
    for n in ast.walk(definition):
        assert not isinstance(n,(ast.Import,ast.ImportFrom,ast.Global,ast.Nonlocal,ast.With,ast.Try))
        if isinstance(n,ast.Call):
            assert (isinstance(n.func,ast.Name) and n.func.id in allowed) or (isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='np' and n.func.attr in ['maximum','sqrt'])
reference={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(curves),'exec'),reference)
exact_nodes=[n for n in ast.parse(probe.read_text(encoding='utf8')).body if isinstance(n,ast.FunctionDef) and n.name=='exact']
assert len(exact_nodes)==1
for n in ast.walk(exact_nodes[0]):
    assert not isinstance(n,(ast.Import,ast.ImportFrom,ast.Global,ast.Nonlocal,ast.With,ast.Try))
    if isinstance(n,ast.Call) and isinstance(n.func,ast.Name):
        assert n.func.id in {'ValueError','wave','equation','phiL','phiR','brentq','bisect','max','min','abs','float','density','all','dict','list'}
namespace=dict(np=np,math=math,brentq=brentq,bisect=bisect,reference=reference)
exec(compile(ast.Module(body=exact_nodes,type_ignores=[]),str(probe),'exec'),namespace)
s=loadmat(file,simplify_cells=True);g=s['gas'];n=s['normal'];report=s['report']
states=[]
for key in ['ql','qr']:
    q=s[key];rho=q[0];u=q[1:4]/rho;internal=q[4]-.5*rho*np.dot(u,u)
    pressure=(g['gamma']-1)*internal;temp=pressure/(rho*g['R']);assert rho>0 and internal>0
    assert g['temperatureRange'][0]<=temp<=g['temperatureRange'][1]
    states.append(dict(rhoKgM3=float(rho),normalVelocityMs=float(u@n),pressurePa=float(pressure),temperatureK=float(temp)))
left,right=states
answer=namespace['exact']((left['rhoKgM3'],left['normalVelocityMs'],left['pressurePa']),
                          (right['rhoKgM3'],right['normalVelocityMs'],right['pressurePa']),g['gamma'])
rho,un,pressure=answer['sample'];temperature=pressure/(rho*g['R'])
recorded=float(re.search(r'T=([-+0-9.eE]+)',report['message']).group(1))
assert abs(temperature-recorded)<1e-9
assert not(g['temperatureRange'][0]<=temperature<=g['temperatureRange'][1]) and not report['candidateAccepted']
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-rejected-interface-pressure-root-and-wave-sampling',case=case,processID=os.getpid(),
         faceID=int(report['globalInteriorFaceID']),left=left,right=right,exactRiemann=answer,
         independentInterfaceTemperatureK=float(temperature),recordedInterfaceTemperatureK=recorded,
         interfaceTemperatureDifferenceK=float(temperature-recorded),originalTemperatureRangeK=g['temperatureRange'].tolist(),
         candidateAccepted=False,fullNSRecomputed=False,physicalSolution=False,newProductionSolvedAreaM2=0,
         inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
