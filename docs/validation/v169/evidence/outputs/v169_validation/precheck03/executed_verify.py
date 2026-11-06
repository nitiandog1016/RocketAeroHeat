"""模块用途：独立验证真实鼻尖、三维混合几何、分段热边界和原文局部测量算子。"""
import argparse
import copy
from pathlib import Path
import sys
import time

import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from naca3344Model import D,R,SHOULDER,END,source_data,wall_temperature,strip_mean,heat_parameter,recovery_temperature
from su2NACALaminarInterface import mesh,validate
from su2ThermalRuntime import read,sha
from nacaWallFlux import geometry,evaluate
from validateIdealAir import save


def hexa_jacobian(q,coordinate):
    """独立标准八节点形函数，不调用生成器的几何/体积算法。"""
    signs=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                    [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
    grad=np.empty((8,3))
    for j in range(3):
        others=[i for i in range(3) if i!=j]
        grad[:,j]=signs[:,j]*np.prod(1+signs[:,others]*coordinate[others],axis=1)/8
    return np.linalg.det(q.T@grad)


def verify(folder,out):
    """模块一：有区分力的入口检查通过后才允许原生三更新。"""
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    (out/'executed_verify.py').write_bytes(Path(__file__).read_bytes())
    p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
    report=read(folder/'geometry_audit.json')
    locks={str(x):sha(x) for x in [folder/'mesh.su2',folder/'geometry.json',folder/'case.cfg',
        folder/'case_manifest.json',folder/'initial_flow.csv',Path(__file__).resolve(),
        ROOT/'scripts/nacaWallFlux.py',ROOT/'scripts/naca3344Model.py']}
    # 所有六面体均是显式角度/法向张量积，正Jacobian下界解析成立，另以Gauss积分复核体积。
    labels=np.array(meta['nodeIndices']);s=np.array(meta['s']);offset=np.array(meta['offsets']);az=np.array(meta['azimuth'])
    lower=[];volume=0.;checkerr=[]
    gauss=np.polynomial.legendre.leggauss(3)
    for c in cells:
        if c[0]!=12:continue
        q=p[c[1:]];j,i,k=labels[c[1]]
        dy=offset[j+1]-offset[j];dr=R+offset[j];da=az[k+1]-az[k]
        if s[i+1]<=SHOULDER+1e-14:
            bound=dy*dr**2*np.sin(s[i]/R)*np.sin((s[i+1]-s[i])/R)*np.sin(da)/8
        else:
            bound=(s[i+1]-s[i])*dy*dr*np.sin(da)/8
        lower.append(bound)
        observed=min(hexa_jacobian(q,np.array([a,b,d],float)) for a in (-1,1) for b in (-1,1) for d in (-1,1))
        checkerr.append(abs(observed/bound-1))
        for a,wa in zip(*gauss):
            for b,wb in zip(*gauss):
                for d,wd in zip(*gauss):volume+=wa*wb*wd*hexa_jacobian(q,np.array([a,b,d]))
    assert min(lower)>0 and max(checkerr)<2e-9
    # 棱柱的三角/层插值，体积分用独立3×3 Gauss（三角形Duffy变换）。
    for c in cells:
        if c[0]!=13:continue
        q=p[c[1:]];bottom=q[:3];delta=q[3:]-bottom
        for a,wa in zip(*gauss):
            u=(a+1)/2
            for b,wb in zip(*gauss):
                v=(b+1)/2*(1-u)
                for d,wd in zip(*gauss):
                    z=(d+1)/2
                    ds=bottom[1]-bottom[0]+z*(delta[1]-delta[0])
                    dt=bottom[2]-bottom[0]+z*(delta[2]-delta[0])
                    dz=delta[0]*(1-u-v)+delta[1]*u+delta[2]*v
                    jac=np.dot(ds,np.cross(dt,dz));assert jac>0
                    volume+=wa*wb*wd*(1-u)*jac/8
    volume_error=abs(volume/report['fluidVolume_m3']-1);assert volume_error<2e-10
    assert np.array_equal(p[0],[0,0,0]) and np.all(p>=np.array([-.127-1e-15,0,0]))
    assert set(np.unique([c[0] for c in cells]))=={12,13}
    # 段Tw及延拓点明确，面积/符号独立解析例；原测量点本身不允许称W/m²。
    dataset,tw,heat=source_data()
    actual,extended=wall_temperature(tw['s_over_D']*D)
    assert np.max(abs(actual-tw['Tw_K']))<1e-12 and not extended.any()
    assert wall_temperature(END)[0]==tw['Tw_K'][-1] and bool(wall_temperature(END)[1])
    mean,area=strip_mean(lambda x:2.,0,SHOULDER)
    assert abs(mean-2)<1e-12 and abs(area/(2*np.pi*R*R)-1)<1e-12
    # 热通量双marker交界受不同Tw影响，逐边积分合计；以唯一节点值覆盖应产生可见偏差。
    geom=geometry(p,cells,markers)
    temp=np.full(len(p),310.);conductivity=np.full(len(p),.03)
    wall_values={row['marker']:row['givenTw_K'] for row in meta['walls']}
    calculated=evaluate(geom,temp,conductivity,wall_values)
    by_node={}
    for tag,value in calculated.items():
        for node,q,a in zip(value['ids'],value['qIntoWall'],value['dualArea']):
            by_node.setdefault(int(node),[]).append((float(q),float(a)))
    shared=[row for row in by_node.values() if len(row)>1]
    assert shared and all(q<0 for row in by_node.values() for q,a in row)
    full=sum(v['powerIntoWall'] for v in calculated.values())
    overwritten=sum(row[-1][0]*sum(a for q,a in row) for row in by_node.values())
    assert abs(full-overwritten)>1
    # 解析已知输入：Re=100000，Nu=100，参数1/sqrt(10)；错号/零驻点/来流标识拒绝。
    edge=dict(temperatureK=200.,densityKgM3=1.,speedMS=100.,viscosityPas=1e-5,conductivityWmK=.03)
    r=heat_parameter(.01,-15000.,350.,edge,300.,dict(edge='synthetic-local-state',recovery='synthetic-analytic'))
    assert abs(float(r['parameter'])-1/np.sqrt(10))<1e-14
    assert float(recovery_temperature(200.,320.,.8,'synthetic'))==296.
    rejected=[]
    for name,call in [('wrong-sign',lambda:heat_parameter(.01,15000.,350.,edge,300.,dict(edge='e',recovery='r'))),
        ('stagnation-0/0',lambda:heat_parameter(0.,-15000.,350.,edge,300.,dict(edge='e',recovery='r'))),
        ('missing-edge',lambda:heat_parameter(.01,-15000.,350.,edge,300.,{})),
        ('unheated-nose-Tr',lambda:recovery_temperature(200.,320.,.8,'unheated-wall-near-nose'))]:
        try:call()
        except AssertionError:rejected.append(name)
        else:raise AssertionError('应拒绝：'+name)
    manifest=read(folder/'case_manifest.json');text=(folder/'case.cfg').read_text()
    runtime=read(ROOT/'tools/runtime/su2'/('8.6.0-'+manifest['requiredRuntimeFlavor'])/'runtime.json')
    validate(folder,manifest,text,runtime)
    for name,wrong in [('axis-enabled',text.replace('AXISYMMETRIC= NO','AXISYMMETRIC= YES')),
        ('wrong-Tw',text.replace('wall_000, ','wall_000, 0 # ',1)),
        ('RANS',text.replace('SOLVER= NAVIER_STOKES','SOLVER= RANS')),
        ('iterations-unplanned',text.replace('ITER= 3\n','ITER= 3000\n'))]:
        try:validate(folder,manifest,wrong,runtime)
        except AssertionError:rejected.append(name)
        else:raise AssertionError('应拒绝：'+name)
    from su2IdealAirInterface import settings,OBSERVED_FLAVOR
    if manifest['requiredRuntimeFlavor']==OBSERVED_FLAVOR:
        tokens=settings(text)['MARKER_RIEMANN'].strip('() ').split(',')
        mach_vector=np.array(tokens[4:7],float)
        assert np.array_equal(mach_vector,[-dataset['freeStream']['Mach'],0,0])
        rejected.append('unit-direction-is-not-Mach-vector')
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    result=dict(passed=True,nodes=len(p),cells=len(cells),hexaJacobianMinimumAnalyticLowerBound=min(lower),
        hexaBoundIndependentMaximumError=max(checkerr),volumeIndependentRelativeError=volume_error,
        actualNoseIncluded=True,sharedMarkerNodes=len(shared),syntheticPowerIntoWallW=full,
        wrongOverwritePowerDifferenceW=overwritten-full,rejected=rejected,thermalPoints=len(tw),heatPoints=len(heat),
        inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started,
        newPhysicalQualification=False,readyForRegisteredNativeSmallCheck=True)
    save(out/'report.json',result);print(result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True)
    parser.add_argument('--output',required=True);args=parser.parse_args()
    verify(Path(args.folder).resolve(),Path(args.output).resolve())
