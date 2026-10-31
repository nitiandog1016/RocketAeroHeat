"""最小可区分的状态/几何负例，防止k能量约定和两种壁距混淆。"""
from pathlib import Path
import sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2SSTInterface import state,triangle_distances
from prepareSU2SSTCase import prepare
from su2ThermalRuntime import validate_sst_manifest

class ContractTests(unittest.TestCase):
    def test_named_energy_and_k(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'r.csv'
            p.write_text('"PointID","x","y","Density","Momentum_x","Momentum_y","Energy","Turb_Kin_Energy","Omega"\n0,0,0,2,6,8,225,50,7\n1,1,0,2,6,8,225,50,7\n')
            _,_,nd,r,u,e=state(p,True);self.assertEqual(nd,2);np.testing.assert_allclose(e,50)
            with self.assertRaises(AssertionError):state(p,False)
            p.write_text(p.read_text().replace('225,50,7','125,60,7'))
            with self.assertRaises(AssertionError):state(p,True)
    def test_triangle_projection_interior_and_outside(self):
        a=np.array([[0.,0,0]]);b=np.array([[1.,0,1]]);c=np.array([[0.,1,0]])
        n=np.array([-1.,0,1])/np.sqrt(2);q=np.array([.25,.25,.25])+n*.02
        self.assertAlmostEqual(triangle_distances(q,a,b,c)[0],.02,places=12)
        self.assertAlmostEqual(triangle_distances([0,2,0],a,b,c)[0],1,places=12)
    def test_rans_without_contract_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises((AssertionError,FileNotFoundError,KeyError)):
                validate_sst_manifest(Path(d),{},'SOLVER= RANS\n')
    def test_missing_provenance_before_mutation(self):
        with self.assertRaises(AssertionError):prepare(ROOT/'outputs/v124_validation/missing')

if __name__=='__main__':unittest.main()
