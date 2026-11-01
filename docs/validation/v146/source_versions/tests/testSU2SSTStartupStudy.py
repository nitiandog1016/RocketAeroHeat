"""启动研究的差异约束与恢复许可反例；不重复既有SST格式测试。"""
import copy,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareSU2SSTStartupStudy import KIND,STUDY,original_source,config,settings,changes,validate_comparison
from compareSU2SSTStartupStudy import RESIDUALS,restore_gates

class StartupStudy(unittest.TestCase):
    def test_single_numerical_change(self):
        source=config(original_source()/'case.cfg')
        for variant in ('control','first_order','restore_second_order'):
            old=dict(source)
            if variant=='restore_second_order':old['MUSCL_FLOW']='NO'
            init=dict(kind=KIND,studyID=STUDY,variant=variant,sameMesh=True,samePhysics=True,
                localIterationStartsAt=0,sourceLocalIteration=99 if variant=='restore_second_order' else 149,
                physicsQualified=False,uninterruptedTrajectoryEquivalent=False,sourceStateHasNativePrimitive=True,
                restartCFLReset=True,numericalChanges=changes(old,variant),samePhysicsAndNumerics=variant=='control',
                flowSpatialOrder=1 if variant=='first_order' else 2)
            target=settings(old,variant);validate_comparison(old,target,init)
            for key,value in [('CFL_NUMBER','0.5'),('CONV_RESIDUAL_MINVAL','-6'),('MUSCL_TURB','YES'),('ITER','200'),('NUM_METHOD_GRAD','GREEN_GAUSS')]:
                bad=dict(target);bad[key]=value
                with self.assertRaises(AssertionError):validate_comparison(old,bad,init)
            bad=copy.deepcopy(init);bad['sourceLocalIteration']=75
            with self.assertRaises(AssertionError):validate_comparison(old,target,bad)

    def test_each_failed_gate_blocks_restoration(self):
        source=dict(kFloorOffWall=100,enthalpyExcess=.4)
        control=dict(interfacePassed=True,last20MeanLogResidual={k:1. for k in RESIDUALS},kFloorOffWall=100,
            omegaFloor=0,enthalpyExcess=.4,heatRelativeL1Change=.2,last25HeatPowerRelativeChange=-.2)
        candidate=copy.deepcopy(control);candidate.update(kFloorOffWall=40,enthalpyExcess=.2,heatRelativeL1Change=.1,last25HeatPowerRelativeChange=.1)
        candidate['last20MeanLogResidual']['rms[RhoE]']=.5
        self.assertTrue(all(restore_gates(source,control,candidate)[0].values()))
        for key,value in [('interfacePassed',False),('kFloorOffWall',51),('omegaFloor',1),('enthalpyExcess',.31),('heatRelativeL1Change',.21),('last25HeatPowerRelativeChange',-.21)]:
            bad=copy.deepcopy(candidate);bad[key]=value
            self.assertFalse(all(restore_gates(source,control,bad)[0].values()),key)
        for key,value in [('rms[RhoE]',.8),('rms[k]',1.4),('rms[Rho]',1.4),('rms[w]',1.4)]:
            bad=copy.deepcopy(candidate);bad['last20MeanLogResidual'][key]=value
            self.assertFalse(all(restore_gates(source,control,bad)[0].values()),key)

if __name__=='__main__':unittest.main()
