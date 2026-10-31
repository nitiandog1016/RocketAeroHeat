"""独立复核内部粗解到细网格初场；采用SciPy插值器对照实现而不重算流场。"""
from pathlib import Path
import argparse
import importlib.util
import hashlib
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('target');args=parser.parse_args()
    target=Path(args.target).resolve();report=json.loads((target/'warm_start_manifest.json').read_text())
    for name,h in report['inputSHA256'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h
    source=Path(report['source']);raw=np.loadtxt(str(source/'restart_flow.csv'),delimiter=',',skiprows=1)
    new=np.loadtxt(str(target/'initial_flow.csv'),delimiter=',',skiprows=1)
    # 通过CSV本身恢复张量网格，与生成器读取SU2网格的索引方式独立。
    x=np.unique(raw[:,1]);y=np.unique(raw[:,2]);order=np.lexsort((raw[:,2],raw[:,1]))
    grid=raw[order,3:].reshape(len(x),len(y),4)
    query=new[:,1:3].copy();query[:,0]=np.clip(query[:,0],x[0],x[-1]);query[:,1]=np.clip(query[:,1],y[0],y[-1])
    expected=RegularGridInterpolator((x,y),grid,bounds_error=True)(query)
    relative=np.max(np.abs(expected-new[:,3:])/np.maximum(1,np.abs(expected)),axis=0)
    assert np.max(relative)<2e-12
    rho=new[:,3];vel=new[:,4:6]/rho[:,None];T=(new[:,6]/rho-(vel**2).sum(axis=1)/2)*.4/287.05287
    assert np.all(rho>0) and np.all(T>0)
    wall=(new[:,2]==0)&(new[:,1]<=0);assert np.max(np.abs(vel[wall]))<1e-9
    # 非均匀网格上的双线性制造场应逐点精确，再检查拒绝真实外推。
    path=ROOT/'scripts/prepareSU2FlatPlateWarmStart.py';spec=importlib.util.spec_from_file_location('warm',str(path))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    xx=np.array([0.,.1,.6,1.]);yy=np.array([0.,.2,.8]);X,Y=np.meshgrid(xx,yy,indexing='ij')
    g=np.stack((1+X+Y+X*Y,2-2*X+3*Y+X*Y),axis=2);points=np.array([[.05,.15],[.42,.52],[1.,.8]])
    actual,_=m.interpolate(xx,yy,g,points);a,b=points.T
    assert np.max(abs(actual-np.column_stack((1+a+b+a*b,2-2*a+3*b+a*b))))<2e-15
    try:m.interpolate(xx,yy,g,np.array([[1.01,.2]]))
    except AssertionError:pass
    else:raise AssertionError('外推未被拒绝')
    answer=dict(independentReference='SciPy RegularGridInterpolator on original CSV Cartesian nodes',
                maximumRelativeDifference=relative.tolist(),positiveState=True,noSlipPreserved=True,
                analyticalBilinearCheck=True,extrapolationRejected=True,physicsQualified=False,
                sourceSHA256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                              (path,Path(__file__).resolve(),target/'initial_flow.csv',target/'warm_start_manifest.json')})
    (target/'initial_audit.json').write_text(json.dumps(answer,indent=2)+'\n',encoding='utf8')
    print(json.dumps(answer['maximumRelativeDifference']))


if __name__=='__main__':main()
