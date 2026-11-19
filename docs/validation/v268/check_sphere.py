"""独立用整域配点解核对等温微球的有限 Pe 热流，不读取分区求解器的速度或接口值。"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre
from scipy.io import loadmat
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import spsolve


def derivative(x, weights):
    difference = x[:, None] - x[None, :]
    np.fill_diagonal(difference, 1)
    d = weights[None, :] / weights[:, None] / difference
    np.fill_diagonal(d, 0)
    np.fill_diagonal(d, -d.sum(axis=1))
    return d


def whole_shell(nr, nm, bath, pe):
    # 使用公开整球 Stokes 解直接构造速度，导热和对流在连续方程中分别离散。
    r = 1 + (bath-1)*(1-np.cos(np.arange(nr)*np.pi/(nr-1)))/2
    mu, w = legendre.leggauss(nm)
    wr = (-1.)**np.arange(nr)
    wr[[0,-1]] *= .5
    dr = derivative(r,wr)
    dm = derivative(mu,(-1.)**np.arange(nm)*np.sqrt(w*(1-mu**2)))
    radial = dr@dr + (2/r)[:,None]*dr
    angular = (1-mu**2)[:,None]*(dm@dm)-2*mu[:,None]*dm
    u_r = (1-1.5/r+.5/r**3)[:,None]*mu
    u_mu = ((1-.75/r-.25/r**3)/r)[:,None]*(1-mu**2)
    op = (kron(radial,eye(nm)) + kron(diags(1/r**2),angular) -
          pe*(diags(u_r.ravel())@kron(dr,eye(nm)) + diags(u_mu.ravel())@kron(eye(nr),dm))).tolil()
    rhs = np.zeros(nr*nm)
    for j in range(nm):
        for i in (j,(nr-1)*nm+j):
            op.rows[i]=[i]
            op.data[i]=[1.]
        rhs[j]=1
    theta=spsolve(op.tocsr(),rhs).reshape(nr,nm)
    flux=-(dr@theta)[0]
    return legendre.legfit(mu,flux,nm-1),np.dot(w,flux)


def read_bands(path):
    with path.open(encoding='utf8',newline='') as f:
        rows=list(csv.DictReader(f))
    return [np.array([float(row[key]) for row in rows]) for key in ('muLow','muHigh','heatIntoWallWm2')]


def verify(folder):
    process=json.loads((folder/'debugger_process.json').read_bytes())
    assert process['status']=='terminal' and process['debuggerExitCode']==0 and not process['survivingProcessIDs']
    assert process['stopReason'] is None and all(process['inputsStable'].values())
    record=json.loads((folder/'result.json').read_bytes())
    cases=record['cases'];c=record['condition'];radius=c['radiusM']
    speed=np.linalg.norm(c['velocityMps'])
    pe=c['densityKgM3']*c['specificHeatJkgK']*speed*radius/c['conductivityWmK']
    scale=c['conductivityWmK']*(record['wallTemperatureK']-c['ambientTemperatureK'])/radius
    poly1,power1=whole_shell(33,12,record['bathRatio'],pe)
    poly,power=whole_shell(49,16,record['bathRatio'],pe)
    points=np.linspace(-1,1,101)
    reference_change=max(abs(legendre.legval(points,poly)-legendre.legval(points,poly1)))
    assert reference_change<1e-7 and abs(power1/power-1)<1e-7
    integral=legendre.legint(poly)
    reference_w=-power*2*np.pi*radius**2*scale
    errors={}
    for name in ('coarse','medium','fine','moved'):
        lo,hi,wall=read_bands(folder/name/'wall_bands.csv')
        averaged=(legendre.legval(hi,integral)-legendre.legval(lo,integral))/(hi-lo)
        errors[name]=float(max(abs(wall/scale+averaged)))
        tabulated=float(np.dot(wall,hi-lo)*2*np.pi*radius**2)
        assert math.isclose(tabulated,cases[name]['interface']['curvedWallPowerW'],rel_tol=1e-12)
    assert errors['fine']<.005 and errors['fine']<errors['medium']<errors['coarse']
    assert abs(cases['fine']['interface']['curvedWallPowerW']/reference_w-1)<.002
    _,_,fine=read_bands(folder/'fine/wall_bands.csv')
    _,_,moved=read_bands(folder/'moved/wall_bands.csv')
    cut_change=float(max(abs(fine-moved))/scale)
    assert cut_change<.005 and abs(cases['moved']['interface']['curvedWallPowerW']/cases['fine']['interface']['curvedWallPowerW']-1)<.002
    assert fine[0]<fine[-1]<0
    # 用工作台返回的原节点、三角形重新计算面积，并逐原面重建两个参考场。
    js=(folder/'fine/workbench/case-data.js').read_text(encoding='utf8').strip()
    assert js.startswith('window.ROCKET_CASE=') and js.endswith(';')
    data=json.loads(js[len('window.ROCKET_CASE='):-1])
    vertices=np.asarray(data['vertices']);faces=np.asarray(data['faces']);tri=vertices[faces]
    area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    centers=tri.mean(axis=1)-np.array(c['centerM'])
    mu=(centers @ (np.array(c['velocityMps'])/speed))/np.linalg.norm(centers,axis=1)
    fields={f['key']:f for f in data['fields']}
    p=np.array(fields['pressure']['values'])*1000
    q=np.array(fields['heatFlux']['values'])*1000
    assert all(v==1 for key in ('pressure','heatFlux') for v in fields[key]['status'])
    delta=c['viscosityPaS']*speed/radius
    assert max(abs((p-c['ambientPressurePa'])/delta+1.5*mu))<1e-6
    lo,hi,_=read_bands(folder/'fine/wall_bands.csv')
    match=np.searchsorted(hi,mu,side='right').clip(0,len(hi)-1)
    np.testing.assert_allclose(q,fine[match],atol=1e-10,rtol=1e-12)
    assert math.isclose(float(np.dot(q,area)),cases['fine']['wallTrianglePowerW'],rel_tol=1e-12)
    state=loadmat(folder/'fine/terminal.mat',simplify_cells=True,variable_names=['ok','thermal','mInner','mOuter','advected'])
    assert state['ok'] and np.allclose(state['mInner'],state['mOuter'],rtol=1e-8,atol=1e-23)
    actual=cases['fine']['interface']
    close=lambda a,b:math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-18)
    assert close(sum(state['mInner'][state['mInner']>0]),actual['interfaceMassOutKgS'])
    assert close(-sum(state['mInner'][state['mInner']<0]),actual['interfaceMassInKgS'])
    assert close(sum(state['advected']),actual['interfaceAdvectedHeatNetOutW'])
    assert close(actual['interfaceAdvectedHeatNetOutW']+actual['interfaceConductedHeatNetOutW'],-actual['curvedWallPowerW'])
    thermal=state['thermal']
    np.testing.assert_allclose(thermal['temperatureK'],c['ambientTemperatureK']+
                               (record['wallTemperatureK']-c['ambientTemperatureK'])*thermal['normalizedTemperature'],atol=1e-12,rtol=0)
    np.testing.assert_allclose(fields['nearWallCellTemperature']['values'],thermal['temperatureK'][1,match],atol=1e-12,rtol=0)
    assert all(fields[k]['values']==[None]*len(area) for k in ('wallSolved','edgeTemperature','recoveryTemperature'))
    max_point_error=float(max(abs(q/scale+legendre.legval(mu,poly))))
    assert max_point_error<.005
    return dict(referenceMethod='unsplit Chebyshev radial / Legendre polar collocation with exact Stokes velocity',
                fineNormalizedBandError=errors['fine'],allNormalizedBandErrors=errors,referenceChange=reference_change,
                cutMoveNormalizedBandChange=cut_change,referenceCurvedPowerW=reference_w,
                computedCurvedPowerW=cases['fine']['interface']['curvedWallPowerW'],originalTrianglePowerW=float(np.dot(q,area)),
                originalFaces=len(area),originalAreaM2=float(area.sum()),physicalSphereAreaM2=4*np.pi*radius**2,
                maximumOriginalFacePointFluxErrorScaled=max_point_error,verifiedNearWallTemperatureFaces=len(area),
                nearWallTemperatureRangeK=[float(np.min(fields['nearWallCellTemperature']['values'])),float(np.max(fields['nearWallCellTemperature']['values']))],
                upstreamBandWm2=float(fine[0]),downstreamBandWm2=float(fine[-1]),
                sourceResultSHA256=hashlib.sha256((folder/'result.json').read_bytes()).hexdigest(),
                experimentalValidation=False,rocketCoverageGainM2=0)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    parser.add_argument('new_report',type=Path)
    args=parser.parse_args();out=verify(args.folder)
    with args.new_report.open('x',encoding='utf8') as f:
        json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(out,ensure_ascii=True))
