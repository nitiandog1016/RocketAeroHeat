"""新诊断身份的负例：不能借单步证据绕过预算、源失败身份或原物理配置。"""
from pathlib import Path
import sys,unittest,copy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2SSTPhaseProbe import KIND,validate_settings

class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.c=dict(ITER='1',OUTPUT_WRT_FREQ='1',OUTPUT_FILES='(SURFACE_PARAVIEW_ASCII, RESTART_ASCII, CSV)',
            VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE)',CONV_RESIDUAL_MINVAL='-8',PRANDTL_TURB='.9')
        self.old=dict(self.c,ITER='100',OUTPUT_WRT_FREQ='50',OUTPUT_FILES='old')
        self.init=dict(kind=KIND,diagnosticOnly=True,sameMesh=True,samePhysicsAndNumerics=True,
            localIterationStartsAt=0,requestedUpdates=1,sourceNativeSucceeded=True,sourcePostAuditPassed=False,physicsQualified=False)
    def test_one_diagnostic_accepted(self):validate_settings(self.c,self.old,self.init)
    def test_multiple_steps_rejected(self):
        for key,value in [('ITER','100'),('OUTPUT_WRT_FREQ','50'),('PRANDTL_TURB','.8'),('CONV_RESIDUAL_MINVAL','-4')]:
            with self.subTest(key=key),self.assertRaises(AssertionError):validate_settings(dict(self.c,**{key:value}),self.old,self.init)
    def test_fake_lineage_rejected(self):
        for key,value in [('kind','same-mesh-laminar-to-SST-initial-guess'),('diagnosticOnly',False),
            ('sourcePostAuditPassed',True),('physicsQualified',True),('localIterationStartsAt',100)]:
            with self.subTest(key=key),self.assertRaises(AssertionError):validate_settings(self.c,self.old,dict(self.init,**{key:value}))

if __name__=='__main__':unittest.main()
