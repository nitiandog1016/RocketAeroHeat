# -*- coding: utf-8 -*-
"""从已绑定公开Clawpack压力曲线及独立Brent采样器取得三维边界参考。"""
from pathlib import Path
import ast,hashlib,math
import numpy as np
from scipy.optimize import brentq,bisect

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';REF=ROOT/'outputs/v82_validation'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

# 优先采用随库的小型原文，便于清理outputs或重新下载源码后离线复现独立核对。
bundled=ROOT/'tests/data/riemann_reference/raw'
source=(bundled if (bundled/'clawpack_exact_euler.py').is_file() else REF)/'clawpack_exact_euler.py'
probe=(bundled if (bundled/'probe_exact_riemann.py').is_file() else REF)/'probe_exact_riemann.py'
assert sha(source)=='2a926d59d4935807b93c425aa617bf64d64858e65231fd556b6ae6202d0ec2ef'
assert sha(probe)=='402ffeee97aed006c6784fa6b46d3276a3ab428897d77d07921308d9708c067a'
allowed={'pospart','sound_speed','beta','integral_curve_1','integral_curve_3','hugoniot_locus_1','hugoniot_locus_3'}
nodes=[n for n in ast.parse(source.read_text(encoding='utf8')).body if isinstance(n,ast.FunctionDef) and n.name in allowed]
assert {n.name for n in nodes}==allowed
reference={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),reference)
nodes=[n for n in ast.parse(probe.read_text(encoding='utf8')).body if isinstance(n,ast.FunctionDef) and n.name=='exact']
assert len(nodes)==1
namespace=dict(np=np,math=math,brentq=brentq,bisect=bisect,reference=reference)
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(probe),'exec'),namespace)
REGIONS=['left-state','left-star','left-rarefaction-fan','right-state','right-star','right-rarefaction-fan']

def exact(p,r,n,gas):
    g=float(gas['gamma']);R=float(gas['R'])
    ul=float(p[1:4]@n);ur=float(r[1:4]@n)
    s=namespace['exact']((float(p[0]),ul,float(p[0]*R*p[4])),(float(r[0]),ur,float(r[0]*R*r[4])),g)
    rho,un,pressure=s['sample'];tangent=p[1:4]-ul*n if s['side']=='left' else r[1:4]-ur*n
    velocity=tangent+un*n;energy=pressure/(g-1)+.5*rho*float(velocity@velocity)
    state=np.r_[rho,velocity,pressure/(rho*R)]
    flux=np.r_[rho*un,rho*velocity*un+pressure*n,(energy+pressure)*un]
    return state,flux,REGIONS.index(s['sampleRegion'])+1

PUBLIC_BINDINGS={str(p.relative_to(ROOT)):sha(p) for p in (source,probe)}
