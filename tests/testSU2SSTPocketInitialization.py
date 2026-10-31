"""检查初猜的解析图延拓、浮点下限及实际字段隔离，不启动原生求解器。"""
from pathlib import Path
import sys,json,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2SSTPocketStudy import reconstruct,floor_nodes,validate_initial
from prepareSU2SSTCheckpoint import validate_source_artifacts
from prepareSU2CheckpointContinuation import read,sha


class PocketTests(unittest.TestCase):
    def test_floor_roundoff_not_physical_threshold_change(self):
        values=np.array([1e-4,1.000000000000001e-4,1.0000001e-4,.1])
        self.assertEqual(floor_nodes(values).tolist(),[0,1])

    def test_analytic_log_extension(self):
        # 中心点0到三个锚的距离均1，故其log值为三锚log的算术平均。
        p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
        ids,w,r=reconstruct(p,[[10,0,1,2,3]],[1,2,3],np.array([1e-4,1,8,64.]))
        self.assertEqual(ids.tolist(),[0]);self.assertAlmostEqual(w[0],8.,places=12)
        self.assertLess(r['linearEquationRelativeError'],1e-12)

    def test_wall_cannot_be_modified(self):
        p=np.eye(4,3)
        with self.assertRaises(AssertionError):reconstruct(p,[[10,0,1,2,3]],[0],np.array([1e-4,1.,2.,3.]))

    def test_actual_only_146_omega_fields_changed(self):
        p=ROOT/'outputs/v128_validation/pocket_repaired_roundoff';init=read(p/'sst_initialization.json')
        validate_initial(p,Path(init['sourceDirectory']),init)
        # 未声明变换的普通原场复用保护仍拒绝该候选，不能绕过来源检查。
        with self.assertRaises(AssertionError):validate_source_artifacts(p,init)


if __name__=='__main__':
    r=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(PocketTests))
    out=ROOT/'outputs/v128_validation/pocket_tests.json';assert not out.exists()
    out.write_text(json.dumps(dict(passed=r.wasSuccessful(),tests=r.testsRun,solverRuns=0,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/prepareSU2SSTPocketStudy.py',
            ROOT/'scripts/prepareSU2SSTCheckpoint.py']}),indent=2)+'\n',encoding='utf-8')
    raise SystemExit(0 if r.wasSuccessful() else 1)
