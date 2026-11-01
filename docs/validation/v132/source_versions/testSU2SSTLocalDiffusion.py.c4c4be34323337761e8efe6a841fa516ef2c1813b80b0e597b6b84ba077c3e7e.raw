"""局部中位对偶、仿射梯度及非正交拆分的独立解析检查。"""
from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from reviewSU2SSTLocalDiffusion import local_geometry,gradients,projections
class LocalDiffusion(unittest.TestCase):
    def setUp(self):
        self.p=np.array([(x,y,z) for x in [-1.,0.,1.] for y in [-1.,0.,1.] for z in [-1.,0.,1.]])
        idx=lambda x,y,z:(x*3+y)*3+z
        cells=[]
        for x in range(2):
            for y in range(2):
                for z in range(2):
                    for t in [[(x,y),(x+1,y),(x+1,y+1)],[(x,y),(x+1,y+1),(x,y+1)]]:
                        cells.append([idx(a,b,z) for a,b in t]+[idx(a,b,z+1) for a,b in t])
        self.g=local_geometry(self.p,[np.array(cells),np.empty((0,4),int)],{13})
    def test_median_geometry_and_affine_field(self):
        g=self.g[13];self.assertAlmostEqual(g['volume'],1.,places=13)
        self.assertLess(g['closure'],1e-13)
        values=(self.p@np.array([1.,2.,3.]))[:,None];gr=gradients(self.g,values)[13]
        np.testing.assert_allclose(gr['wls'][0],[1,2,3],rtol=0,atol=2e-13)
        np.testing.assert_allclose(gr['gg'][0],[1,2,3],rtol=0,atol=2e-13)
    def test_constant_and_quadratic_flux(self):
        row=self.g[13];total=0.
        for j,e,n in zip(row['ids'],row['edges'],row['normals']):
            orth,corr=projections(e,n,np.array([0.]),np.zeros((1,3)))
            self.assertEqual(float(orth+corr),0.)
            orth,corr=projections(e,n,np.array([self.p[j]@self.p[j]]),self.p[j][None,:])
            total+=float(orth+corr)
        self.assertAlmostEqual(total/row['volume'],6.,places=12)
    def test_nonorthogonal_correction_can_decrease_minimum(self):
        a,b=projections(np.array([1.,0,0]),np.array([1.,1,0]),np.array([1.]),np.array([[0.,-10,0]]))
        self.assertEqual(float(a),1.);self.assertEqual(float(b),-10.)
        self.assertLess(float(a+b),0)
if __name__=='__main__':unittest.main()
