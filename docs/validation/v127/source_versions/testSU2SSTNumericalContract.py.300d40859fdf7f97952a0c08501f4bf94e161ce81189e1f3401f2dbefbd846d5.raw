"""数值方案声明不能绕过同物理约束或旧同配置重启动保护。"""
from pathlib import Path
import sys,copy,json,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import read,config,sha
from prepareSU2SSTNumericalStudy import validate_comparison
from prepareSU2SSTCheckpoint import validate_source


class NumericalContract(unittest.TestCase):
    def setUp(self):
        self.source=config(Path(read(ROOT/'outputs/v126_validation/main01/engine_folder.json')['folder'])/'case.cfg')
        self.folder=ROOT/'outputs/v127_validation/cap05'
        self.target=config(self.folder/'case.cfg');self.init=read(self.folder/'sst_initialization.json')

    def test_two_actual_preparations(self):
        for name in ['control','cap05']:
            folder=ROOT/'outputs/v127_validation'/name
            validate_comparison(self.source,config(folder/'case.cfg'),read(folder/'sst_initialization.json'))
        self.assertEqual(sha(ROOT/'outputs/v127_validation/control/initial_flow.csv'),sha(self.folder/'initial_flow.csv'))

    def test_physics_and_unlisted_numerics_rejected(self):
        for key,value in [('PRANDTL_TURB','.8'),('CFL_NUMBER','.05'),('MACH_NUMBER','6'),('LINEAR_SOLVER_ITER','50'),('CFL_ADAPT_PARAM','(0.5, 1.1, 0.05, 1.0)')]:
            with self.subTest(key=key),self.assertRaises(AssertionError):
                validate_comparison(self.source,dict(self.target,**{key:value}),self.init)

    def test_wrong_label_or_undeclared_change_rejected(self):
        for key,value in [('variant','control'),('samePhysicsAndNumerics',True),('numericalChanges',{}),('sourceLocalIteration',48)]:
            init=copy.deepcopy(self.init);init[key]=value
            with self.subTest(key=key),self.assertRaises(AssertionError):validate_comparison(self.source,self.target,init)

    def test_old_continuation_still_rejects_changed_cfl(self):
        init=copy.deepcopy(self.init);init['kind']='same-mesh-SST-checkpoint-research';init['samePhysicsAndNumerics']=True
        with self.assertRaises(AssertionError):validate_source(self.folder,self.target,init)


if __name__=='__main__':
    result=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(NumericalContract))
    p=ROOT/'outputs/v127_validation/numerical_contract_tests.json';assert not p.exists()
    p.write_text(json.dumps(dict(passed=result.wasSuccessful(),tests=result.testsRun,sourceSHA256=sha(Path(__file__).resolve())),indent=2)+'\n',encoding='utf-8')
    raise SystemExit(0 if result.wasSuccessful() else 1)
