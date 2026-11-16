"""模块职责：检验新帧身份、损坏拒绝与逐节点物理能量变换，防止FAS修正冒充残差。"""
from pathlib import Path
import struct,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from fasDevelopmentFieldsV247 import state,residual


class Frames(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=str(ROOT/'outputs/v247_validation'));self.p=Path(self.tmp.name)/'frame.bin'
        self.a=np.zeros((2,34));self.a[:,0]=[1,0];self.a[:,4:6]=1;self.a[:,27:33]=123
    def tearDown(self):self.tmp.cleanup()
    def write(self,a,magic=b'RAH247LE'):
        self.p.write_bytes(struct.pack('<8sQQ',magic,len(a),a.shape[1])+a.astype('<f8').tobytes())
    def test_valid_correction_is_not_residual(self):
        self.write(self.a);self.assertTrue(np.array_equal(state(self.p,2),self.a))
        with self.assertRaises(ValueError):residual(self.p,2)
    def test_legacy_identity_rejected(self):
        self.write(self.a,b'RAH235LE')
        with self.assertRaises(ValueError):state(self.p,2)
    def test_truncation_rejected(self):
        self.write(self.a);self.p.write_bytes(self.p.read_bytes()[:-1])
        with self.assertRaises(ValueError):state(self.p,2)
    def test_duplicate_global_ids_rejected(self):
        self.a[:,0]=0;self.write(self.a)
        with self.assertRaises(ValueError):state(self.p,2)
    def test_initial_unknown_required(self):
        self.write(self.a)
        with self.assertRaises(ValueError):state(self.p,2,True)
        self.a[:,27:33]=np.nan;self.write(self.a);self.assertTrue(np.isnan(state(self.p,2,True)[:,27:33]).all())
    def test_negative_density_rejected(self):
        self.a[0,5]=-1;self.write(self.a)
        with self.assertRaises(ValueError):state(self.p,2)
    def test_energy_converted_before_norm(self):
        a=np.array([[0,2,0,0,0,600003,4],[1,-2,0,0,0,-600005,-4]],float);self.write(a,b'RAH247RS')
        raw,r=residual(self.p,2);self.assertTrue(np.array_equal(r[:,4],[3,-5]));self.assertEqual(float(np.mean(r[:,4]**2)),17)
    def test_nonfinite_residual_rejected(self):
        a=np.zeros((2,7));a[:,0]=[0,1];a[0,4]=np.nan;self.write(a,b'RAH247RS')
        with self.assertRaises(ValueError):residual(self.p,2)


if __name__=='__main__':unittest.main()
