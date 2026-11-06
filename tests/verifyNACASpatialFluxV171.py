"""模块用途：以平面已知功率/分片解析值和既存曲面网格检验空间Tw后处理。"""
import argparse
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from nacaSpatialWallFlux import geometry,evaluate,segment_integrals,read_temperature_table,temperature_at,verify_applied
from analyzeNACASpatialHeat import column_indices,column_edge
from su2NACALaminarInterface import mesh
from su2ThermalRuntime import sha,read
from validateIdealAir import save


def planar(shift=0.,height=1.):
    p=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                [shift,0,height],[1+shift,0,height],[1+shift,1,height],[shift,1,height]],float)
    return p,[[12,0,1,2,3,4,5,6,7]],{'wall':[[9,0,3,2,1]]}


def verify(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    sources=[Path(__file__).resolve(),ROOT/'scripts/nacaSpatialWallFlux.py',ROOT/'scripts/analyzeNACASpatialHeat.py',
             ROOT/'scripts/nacaWallFlux.py',ROOT/'scripts/nacaBoundaryLayerEdge.py']
    locks={str(p):sha(p) for p in sources}
    for p in sources:(out/p.name).write_bytes(p.read_bytes())
    p,c,m=planar();g=geometry(p,c,m,p[:,0]);tw=330+2*p[g['ids'],0]
    t=np.full(len(p),300.);k=np.full(len(p),.02)
    desired=np.array([10.,30.,30.,10.]);t[4:]=tw+desired/.02
    flux=evaluate(g,t,k,tw);assert np.max(abs(flux['qIntoWall']-desired))<1e-12
    assert abs(flux['powerIntoWall']-20)<1e-12
    intervals=np.array([[0.,1.],[.2,.8],[0,.37],[.37,1.]])
    result=segment_integrals(g,flux,intervals,2*np.pi)
    expected=np.array([20.,12.,3.7,16.3])
    error=max(abs(r['powerIntoWallSectorW']-v) for r,v in zip(result,expected))
    area_error=max(abs(r['discreteGeometricAreaSectorM2']-(b-a)) for r,(a,b) in zip(result,intervals))
    assert error<1e-12 and area_error<1e-13
    scaled=segment_integrals(g,flux,[[0,1]],np.pi/2)[0]
    assert abs(scaled['symmetryScaledPowerIntoWallFullW']-80)<1e-12
    # 改变符号仍守恒，不能用abs(q)通过加热实验或积分。
    t[4:]=tw+np.array([-10.,30.,30.,-10.])/.02
    signed=evaluate(g,t,k,tw);assert abs(signed['powerIntoWall']-10)<1e-12
    # 倾斜近乎切向的真实棱边触发原生5%下限，精确投影不是可替代公式。
    skew,cc,mm=planar(100.,.001);gg=geometry(skew,cc,mm,skew[:,0])
    assert gg['distanceFloorTriggered'].all()
    expected_distance=.05*np.sqrt(100.**2+.001**2)
    assert np.max(abs(gg['distance']-expected_distance))<1e-12
    ft=evaluate(gg,np.r_[np.full(4,300.),np.full(4,310.)],np.full(8,.02),np.full(4,300.))
    assert np.max(abs(ft['qIntoWall']-.2/expected_distance))<1e-13
    # 原生表语义：两个端点/中点精确，越界不得np.interp静默钳位。
    table=out/'table.csv';table.write_text('RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K\n0,330\n1,332\n',encoding='ascii')
    table_data=read_temperature_table(table);assert np.array_equal(temperature_at([0,.5,1],table_data),[330,331,332])
    applied=out/'applied.csv'
    applied.write_text('globalNodeID,X_m,Y_m,Z_m,givenTemperature_K\n'+''.join('%d,%.17g,%.17g,%.17g,%.17g\n'%(i,*p[i],tw[i]) for i in (3,1,0,2)),encoding='ascii')
    assert verify_applied(applied,p,g,tw)['maximumGivenTwDifferenceK']==0
    rejected=[]
    for name,call in [('profile-extrapolation',lambda:temperature_at([-1e-15],table_data)),
                      ('reversed-segment',lambda:segment_integrals(g,flux,[[.7,.3]])),
                      ('outside-wall',lambda:segment_integrals(g,flux,[[-1e-4,.5]])),
                      ('wrong-Tw-shape',lambda:evaluate(g,t,k,np.full(8,330.))),
                      ('mismatched-applied',lambda:verify_applied(applied,p,g,tw+1))]:
        try:call()
        except AssertionError:rejected.append(name)
        else:raise AssertionError('Expected rejection '+name)
    # 复用唯一实际几何，只检查新子面裁切/功率分解，绝不生成第二网格或CFD。
    folder=ROOT/'outputs/v170_validation/mesh02'
    for name in ('mesh.su2','geometry.json','mesh_arrays.npz'):
        locks[str(folder/name)]=sha(folder/name)
    pp,cc,mm=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
    labels,ix=column_indices(meta,pp,folder/'mesh_arrays.npz');ss=np.array(meta['s'])[labels[:,1]]
    curved=geometry(pp,cc,mm,ss);ids=curved['ids'];nn=curved['normalNeighbor']
    assert not curved['distanceFloorTriggered'].any()
    temp=np.full(len(pp),335.);cond=np.full(len(pp),.03)
    target=-1000-8000*ss[ids];temp[nn]=335+target*curved['distance']/.03
    fcurved=evaluate(curved,temp,cond,np.full(len(ids),335.))
    cuts=np.array([0,.0123,.037,.078,.101,max(meta['s'])])
    pieces=segment_integrals(curved,fcurved,np.column_stack((cuts[:-1],cuts[1:])))
    whole=segment_integrals(curved,fcurved,[[cuts[0],cuts[-1]]])[0]
    sum_error=abs(sum(row['powerIntoWallSectorW'] for row in pieces)/fcurved['powerIntoWall']-1)
    assert sum_error<2e-12 and abs(whole['powerIntoWallSectorW']/fcurved['powerIntoWall']-1)<2e-12
    assert abs(sum(row['discreteGeometricAreaSectorM2'] for row in pieces)/curved['discreteGeometricArea']-1)<2e-12
    # 法向适配只检验原编号/形状/拒绝，不用合成字段声称实际外缘已接受。
    density=np.ones(len(pp));values=dict(Density=density,Momentum_x=np.zeros(len(pp)),Momentum_y=np.zeros(len(pp)),
        Momentum_z=np.zeros(len(pp)),Pressure=np.full(len(pp),50000.),Temperature=np.full(len(pp),300.),
        Laminar_Viscosity=np.full(len(pp),1.8e-5))
    thermo=dict(h=np.full(len(pp),300000.),cp=np.full(len(pp),1000.))
    i=20;az=12
    mapped=column_edge(meta,pp,values,thermo,ix,meta['s'][i],az,'synthetic-no-flow',state_kind='synthetic-test')
    assert mapped['originalNodeIDs']==ix[:,i,az].tolist() and not mapped['accepted']
    missing=column_edge(meta,pp,values,thermo,ix,.5*(meta['s'][i]+meta['s'][i+1]),az,'synthetic-no-flow',state_kind='synthetic-test')
    assert missing['status']=='measurement-center-not-an-original-column'
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    report=dict(passed=True,planarIndependentPowerAbsoluteErrorW=float(error),planarIndependentAreaErrorM2=float(area_error),
        triggeredFloorNodes=int(gg['distanceFloorTriggered'].sum()),curvedWallNodes=len(ids),
        curvedPartitionPowerRelativeError=float(sum_error),curvedPatchToNativePowerScaledError=fcurved['patchPowerPartitionScaledError'],
        curvedNativeProjectedAreaM2=curved['nativeProjectedArea'],curvedGeometricAreaM2=curved['discreteGeometricArea'],
        originalColumnMapping=True,missingColumnRejected=True,rejected=rejected,inputSHA256=locks,inputsStable=stable,
        newCFDTasks=0,physicalQualification=False,elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report);save(out/'planar_clipped.json',result);save(out/'curved_clipped.json',pieces)
    print({k:v for k,v in report.items() if k not in ('inputSHA256','inputsStable')})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    verify(Path(parser.parse_args().output).resolve())
