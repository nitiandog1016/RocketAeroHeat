"""模块职责：独立小几何与解析特征方向测试，防止分区法向和输入错误被误当物理结果。"""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from phaseEulerRegionalVolume import topology,extract_region,interface_characteristics
from phaseAwareViscousAir import PhaseAwareViscousAir


def main():
    # 两个体积1/6的四面体共享z=0面。顶点次序刻意相反，法向须由owner推断。
    xyz=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]])
    cells=np.array([[0,1,2,3],[0,2,1,4]])
    wall=np.array([[0,1,3],[1,2,3],[2,0,3],[0,2,4],[2,1,4],[1,0,4]])
    top=topology(xyz,cells,{'wall':wall});up=extract_region(top,np.array([True,False]));down=extract_region(top,np.array([False,True]))
    assert up['summary']['components']==down['summary']['components']==1
    assert abs(up['summary']['volumeM3']-1/6)<1e-14 and abs(down['summary']['volumeM3']-1/6)<1e-14
    assert len(up['interface'])==1 and np.array_equal(up['interface'],down['interface'])
    i=up['interface'][0];k=np.flatnonzero(up['boundary']==i)[0];normal=up['normal'][[k]]
    assert np.array_equal(normal,np.array([[0.,0.,-.5]]))
    air=PhaseAwareViscousAir();s=air.evaluate(10000.,300.);rho=float(s['densityKgm3']);a=float(s['soundSpeedMps']);e=float(s['internalEnergyJkg'])
    cases=[]
    for mach,expected in [(-2.,5),(-.5,4),(0.,1),(.5,1),(2.,0)]:
        velocity=mach*a*normal[0]/np.linalg.norm(normal[0]);u=np.r_[rho,rho*velocity,rho*(e+300000+.5*(velocity@velocity))]
        q=interface_characteristics(np.tile(u,(5,1)),top['faces'][[i]],normal)
        assert q['incoming'][0]==expected and abs(q['normalMach'][0]-mach)<1e-13
        assert bool(q['sonic'][0])==(mach==0)
        cases.append(dict(normalMach=mach,incoming=expected))
    failures=[]
    def reject(label,operation):
        try:operation()
        except ValueError:failures.append(label)
        else:raise AssertionError('Expected rejection: '+label)
    reject('missing physical boundary',lambda:topology(xyz,cells,{'wall':wall[:-1]}))
    reject('duplicate physical boundary',lambda:topology(xyz,cells,{'wall':np.r_[wall,wall[:1]]}))
    reject('three owners',lambda:topology(xyz,np.r_[cells,cells[:1]],{'wall':wall}))
    U=np.tile(u,(5,1));U[:,0]=-1
    reject('negative density',lambda:interface_characteristics(U,top['faces'][[i]],normal))
    U=np.tile(u,(5,1));U[:,4]=np.nan
    reject('nonfinite energy',lambda:interface_characteristics(U,top['faces'][[i]],normal))
    U=np.tile(u,(5,1));U[:,1:4]=0;U[:,4]=rho*(float(air.H[0]-air.R*air.T[0])-1+300000)
    reject('below actual caloric domain',lambda:interface_characteristics(U,top['faces'][[i]],normal))
    print(json.dumps(dict(passed=True,analyticBipyramidVolume=1/3,characteristicCases=cases,rejections=failures)))


if __name__=='__main__':main()
