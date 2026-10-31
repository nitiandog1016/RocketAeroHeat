"""候选运行时显式输入/拒绝条件，避免研究选项进入官方或空间壁温路径。"""
from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2ThermalRuntime import validate_frozen_manifest
class FrozenRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.cfg=(ROOT/'docs/validation/v116/cases/control/case.cfg').read_text(encoding='utf-8')
        self.data=dict(frozenEnergyMode='diagnostic',thermalBoundary=dict(mode='adiabatic'))
    def test_explicit_modes(self):
        for mode in ['diagnostic','apply']:
            self.data['frozenEnergyMode']=mode
            self.assertEqual(validate_frozen_manifest(self.data,self.cfg),mode)
    def test_reject_scope(self):
        for old,new in [('SOLVER= NAVIER_STOKES','SOLVER= RANS'),('MGLEVEL= 0','MGLEVEL= 1'),
                        ('REF_DIMENSIONALIZATION= DIMENSIONAL','REF_DIMENSIONALIZATION= FREESTREAM_VEL_EQ_ONE'),
                        ('MARKER_HEATFLUX= (wall, 0)','MARKER_HEATFLUX= (wall, 100)')]:
            with self.assertRaises(AssertionError):validate_frozen_manifest(self.data,self.cfg.replace(old,new))
        for extra in ['TIME_DOMAIN= YES','MARKER_PYTHON_CUSTOM= (wall)','CONFIG_LIST= (a,b)','SOLVER= NAVIER_STOKES']:
            with self.assertRaises(AssertionError):validate_frozen_manifest(self.data,self.cfg+'\n'+extra)
    def test_reject_unknown(self):
        self.data['frozenEnergyMode']='automatic'
        with self.assertRaises(AssertionError):validate_frozen_manifest(self.data,self.cfg)
    def test_isothermal(self):
        self.data['thermalBoundary']['mode']='isothermal'
        cfg=self.cfg.replace('MARKER_HEATFLUX= (wall, 0)','MARKER_ISOTHERMAL= (wall, 300)')
        self.assertEqual(validate_frozen_manifest(self.data,cfg),'diagnostic')
        with self.assertRaises(AssertionError):validate_frozen_manifest(self.data,cfg.replace('(wall, 300)','(wall, -1)'))
if __name__=='__main__':unittest.main()
