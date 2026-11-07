"""独立测站插值真值和拒绝检查；不调用CFD、不依赖实验拟合通过。"""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from compareEilmerRun14Control import interpolate_station,read_wall,tangential_wall_shear


def test(case,out):
    started=time.monotonic();dest=Path(out);assert not dest.exists()
    # 独立旋转例：切向4Pa、法向100Pa，不允许后者改变u_tau=2m/s。
    t=np.array([[.6,.8],[.6,.8],[1.,0.]])
    n=np.array([[-.8,.6],[-.8,.6],[0.,1.]])
    f=4*t+100*n
    stress,yp=tangential_wall_shear(f,np.array([1.,-1.,1.]),t,np.ones(3),np.full(3,.01),np.full(3,.005))
    assert np.allclose(stress,[4.,-4.,4.],rtol=0,atol=2e-14)
    assert np.allclose(yp,1.,rtol=0,atol=3e-15)
    x=np.array([.01,.04,.08,.10,.11,.14,.18,.219]);block=np.array([0,0,0,0,1,1,1,1])
    y=np.where(block==0,100*x-5,-200*x+30)
    positive=0
    for xi in [.01,.025,.08,.09,.11,.125,.17,.2,.219]:
        expected=100*xi-5 if xi<.1017 else -200*xi+30
        got,why=interpolate_station(x,y,block,xi);assert abs(got-expected)<2e-14;positive+=1
    for xi in [0.,.101,.1017,.105,.23]:
        got,why=interpolate_station(x,y,block,xi);assert got is None
    damaged=y.copy();damaged[1]=np.nan
    assert interpolate_station(x,damaged,block,.02)[0] is None
    assert interpolate_station(x,damaged,block,.04)[0] is None
    rejected=[]
    for name,args in [('nan_coordinate',(x,y,block,np.nan)),('unordered',(x[::-1],y[::-1],block[::-1],.15)),
                      ('duplicate',(np.array([.01,.01,.08,.10,.11,.14,.18,.219]),y,block,.02))]:
        try:interpolate_station(*args)
        except ValueError:rejected.append(name)
        else:raise AssertionError(name)
    # 原生第一份完整载荷仅只读；独立解析整周壁面积与固定热流功率。
    a=read_wall(case,0);R=.0325;dx=.22-.1017;rf=R+dx*np.tan(np.pi/6)
    area=2*np.pi*R*.1017+np.pi*(R+rf)*dx/np.cos(np.pi/6)
    assert abs(a[:,3].sum()/area-1)<2e-13
    assert abs(np.dot(np.full(320,100.),a[:,3])/(100*area)-1)<2e-13
    assert np.all(a[:,7]>0) and np.all(a[:,6]==295.2)
    report=dict(passed=True,linearTruthPoints=positive,unknownAndOutOfSupport=7,rejected=rejected,tangentialProjectionTruthPoints=3,
                nativeWallFaces=320,fullRingAreaM2=area,referenceDataRecomputed=False,
                nativeRuns=0,physicsQualified=False,elapsedSeconds=time.monotonic()-started)
    dest.write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(test(a.case,a.out)))
