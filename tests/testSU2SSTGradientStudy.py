"""真实配置的有限差异约束：错误梯度、重构、下限和源迭代必须拒绝。"""
import copy,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareSU2SSTGradientStudy import SOURCE,KIND,STUDY,config,settings,changes,validate_comparison

class GradientStudy(unittest.TestCase):
    def test_allowed_and_rejected_changes(self):
        source=config(SOURCE/'case.cfg')
        for variant in ('wls','gg'):
            init=dict(kind=KIND,studyID=STUDY,variant=variant,sameMesh=True,samePhysics=True,
                localIterationStartsAt=0,sourceLocalIteration=99,physicsQualified=False,
                uninterruptedTrajectoryEquivalent=False,sourceStateHasNativePrimitive=True,
                numericalChanges=changes(variant),samePhysicsAndNumerics=variant=='wls')
            target=settings(source,variant);validate_comparison(source,target,init)
            for key,bad in [('NUM_METHOD_GRAD_RECON','GREEN_GAUSS'),('ITER','100'),
                            ('CONV_RESIDUAL_MINVAL','-6'),('CFL_NUMBER','0.5')]:
                other=dict(target);other[key]=bad
                with self.assertRaises(AssertionError):validate_comparison(source,other,init)
            other=copy.deepcopy(init);other['sourceLocalIteration']=149
            with self.assertRaises(AssertionError):validate_comparison(source,target,other)

if __name__=='__main__':unittest.main()
