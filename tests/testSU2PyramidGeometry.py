"""金字塔参考映射、体积和反向拒绝；数值差分独立检查解析Jacobian。"""
from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from su2PyramidGeometry import pyramid_geometry

def mapped(p,u):
    r,s,t=u
    return (1-t)*((1-r)*(1-s)*p[0]+r*(1-s)*p[1]+r*s*p[2]+(1-r)*s*p[3])+t*p[4]

class PyramidGeometry(unittest.TestCase):
    def test_regular_and_warped(self):
        for p in [np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],[.5,.5,1]]),
                  np.array([[0.,0,0],[1,0,.1],[1,1,-.05],[0,1,0],[.3,.6,1.2]])]:
            volume,report=pyramid_geometry(p,np.arange(5)[None,:]);total=0.;e=1e-6
            g=(1+np.array([-1,1])/np.sqrt(3))/2
            for r in g:
                for s in g:
                    for t in g:
                        u=np.array([r,s,t]);J=np.column_stack([(mapped(p,u+e*np.eye(3)[i])-mapped(p,u-e*np.eye(3)[i]))/(2*e) for i in range(3)])
                        det=np.linalg.det(J);self.assertGreater(det,0)
                        self.assertGreaterEqual(det/(1-t)**2+1e-9,report['minimumScaledJacobian']);total+=det/8
            self.assertAlmostEqual(volume[0],total,places=9)
    def test_inversion_rejected(self):
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],[.5,.5,-1]])
        with self.assertRaises(AssertionError):pyramid_geometry(p,np.arange(5)[None,:])
    def test_empty_supported(self):
        v,r=pyramid_geometry(np.empty((0,3)),np.empty((0,5),int));self.assertEqual(len(v),0);self.assertIsNone(r['minimumScaledJacobian'])

if __name__=='__main__':unittest.main()
