"""扩散—耗散初猜与独立连续解析极限的网格收敛核对，不运行原生求解器。"""
from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2SSTColumnStudy import column_solution


class ColumnTests(unittest.TestCase):
    def test_independent_continuous_solution_convergence(self):
        errors=[];a=.2;nu=.01;beta=.075;m=(1+np.sqrt(1+24*.09/beta))/2
        for n in (64,128):
            y=np.linspace(0,1,n+1);exactw=6*nu/(beta*(y+a)**2)
            exactk=((y+a)**m-a**(2*m-1)*(y+a)**(1-m))/((1+a)**m-a**(2*m-1)*(1+a)**(1-m))
            k,w,ew,ek,it=column_solution(y,nu,0,1,exactw[0],exactw[-1])
            errors.append([max(abs(w/exactw-1)),max(abs(k-exactk))])
            self.assertLess(ew,1e-10);self.assertLess(ek,1e-10);self.assertTrue(np.all(k[1:]>0))
            self.assertEqual(w[0],exactw[0]);self.assertEqual(w[-1],exactw[-1]);self.assertEqual(k[0],0);self.assertEqual(k[-1],1)
        self.assertTrue(np.all(np.array(errors[0])/errors[1]>3.5))
        self.assertTrue(max(errors[1])<.002)

    def test_nonuniform_large_wall_ratio_positive(self):
        y=np.r_[0,np.cumsum(2e-6*1.28471**np.arange(12))]
        k,w,ew,ek,it=column_solution(y,.012,1e-10,500,2.4e12,2e5)
        self.assertTrue(np.all(k>0) and np.all(w>0));self.assertLess(ew,1e-10);self.assertLess(ek,1e-10)

    def test_invalid_geometry_and_boundary_rejected(self):
        for y in ([0,0,1],[0,1],[1,2,3]):
            with self.assertRaises(AssertionError):column_solution(y,.01,0,1,100,1)
        with self.assertRaises(AssertionError):column_solution([0,.5,1],.01,0,1,-1,1)

    def test_native_floor_complementarity(self):
        y=np.linspace(0,.0001,13)
        k,w,ew,ek,it=column_solution(y,.01,1e-10,1e-10,1e12,1e7,1e-10)
        self.assertTrue(np.all(k==1e-10));self.assertLess(ek,1e-10)
        # 两端同为原下限时，耗散缺陷为正而非伪称零方程残差。
        self.assertTrue(np.all(.09*w[1:-1]*k[1:-1]>0))


if __name__=='__main__':unittest.main()
