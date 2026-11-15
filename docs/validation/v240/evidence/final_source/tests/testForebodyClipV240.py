"""模块职责：独立凸包体积和线性场验证裁切，检查真实相邻单元共形及平面边界。"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from scipy.spatial import ConvexHull
from clipEulerForebody import clip_tetrahedra,interpolate_saved_field,write_su2
from phaseEulerRegionalVolume import topology,extract_region,read_tetra_mesh
from runSU2RocketPhaseV222 import save


def boundary(cells):
    faces=np.concatenate([cells[:,j] for j in ([0,1,2],[0,3,1],[1,3,2],[2,3,0])])
    _,first,count=np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_counts=True)
    return faces[first[count==1]]


def oracle(xyz,cut):
    """独立路径：逐原边求交，Qhull算凸包体积，不调用生产三角化或来源权重。"""
    points=[v for v in xyz if v[0]>=cut]
    for i in range(4):
        for j in range(i):
            a,b=xyz[i],xyz[j]
            if (a[0]-cut)*(b[0]-cut)<0:
                points.append(((b[0]-cut)*a+(cut-a[0])*b)/(b[0]-a[0]))
    if len(points)<4:return 0.
    return ConvexHull(np.array(points)).volume


def run(folder):
    folder.mkdir(exist_ok=False)
    xyz=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
    cells=np.array([[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]],np.int64)
    cases=[]
    # 不同顶点顺序和与平面重合的顶点，覆盖1/2/3个原顶点保留的拓扑。
    for name,x,c,cut in [('cube',xyz,cells,.37),('one',np.array([[0,0,0],[2,0,0],[0,1,0],[0,0,1.]]),np.array([[0,1,2,3]]),.5),
      ('two',np.array([[-1,0,0],[2,0,0],[1,1,0],[-2,0,1.]]),np.array([[0,1,2,3]]),0.),
      ('three',np.array([[0,0,0],[2,0,0],[2,1,0],[2,0,1.]]),np.array([[0,1,2,3]]),.5),
      ('on_vertex',np.array([[-1,0,0],[2,0,0],[0,1,0],[0,0,1.]]),np.array([[0,1,2,3]]),0.)]:
        mesh=clip_tetrahedra(x,c,{'farfield':boundary(c)},cut)
        topo=topology(mesh['xyz'],mesh['cells'],mesh['markers']);closed=extract_region(topo,np.ones(len(mesh['cells']),bool))
        expected=sum(oracle(x[row],cut) for row in c)
        error=abs(sum(topo['volume'])-expected)/expected;assert error<2e-14,(name,error)
        if name=='cube':assert abs(expected-.63)<1e-14
        field=np.column_stack([2+x[:,0]-3*x[:,1]+5*x[:,2],x])
        actual=interpolate_saved_field(mesh,field);q=mesh['xyz'];truth=np.column_stack([2+q[:,0]-3*q[:,1]+5*q[:,2],q])
        linear=float(np.max(abs(actual-truth)));assert linear<4e-15,(name,linear)
        normal=topo['normal'][topo['markers']['outlet']];unit=normal/np.linalg.norm(normal,axis=1)[:,None]
        assert np.max(abs(unit-[-1.,0.,0.]))<1e-14
        path=folder/(name+'.su2');write_su2(path,mesh);xx,cc,mm=read_tetra_mesh(path)
        assert np.array_equal(xx,mesh['xyz']) and np.array_equal(cc,mesh['cells'])
        # 单元编号/局部顺序改变仍共享相同多边形对角线；闭合不是仅总面积抵消。
        reverse=clip_tetrahedra(x,c[::-1,::-1],{'farfield':boundary(c)},cut)
        rt=topology(reverse['xyz'],reverse['cells'],reverse['markers'])
        assert abs(sum(rt['volume'])-expected)<2e-14
        cases.append(dict(name=name,volumeM3=float(expected),relativeVolumeError=float(error),linearError=linear,
          cells=len(mesh['cells']),nodes=len(mesh['xyz']),components=closed['summary']['components']))
    rejection=[]
    for name,cut in [('outside',2.),('nan',float('nan'))]:
        try:clip_tetrahedra(xyz,cells,{'farfield':boundary(cells)},cut)
        except ValueError as e:rejection.append(dict(name=name,error=str(e)))
        else:raise AssertionError(name)
    save(folder/'report.json',dict(passed=True,cases=cases,rejections=rejection,physicalQualification=False))
    print('PASS: independent clipped-volume/linear-field/shared-face controls')


if __name__=='__main__':run(Path(sys.argv[1]))
