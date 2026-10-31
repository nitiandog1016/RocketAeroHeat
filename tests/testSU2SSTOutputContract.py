"""验证新格式的真实连续阶段、重启动与拒绝反例，不重新运行求解器。"""
from pathlib import Path
import sys,json,unittest,tempfile,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2SSTInterface import state,phase_columns
from prepareSU2CheckpointContinuation import read,sha,config
from prepareSU2SSTCheckpoint import validate_source

class OutputTests(unittest.TestCase):
    def test_real_saved_stage_and_restart(self):
        errors=[];changes=[]
        for folder,previous in [('format01','restart_flow_000001.csv'),('restart02','initial_flow.csv')]:
            p=ROOT/'outputs/v126_validation'/folder
            # restart02仅保存末局部1；比较前驱需该目录局部0，而原生未输出0。
            if folder=='restart02':
                r=read(p/'engine_process.json');self.assertEqual(r['exitCode'],0)
                self.assertEqual(sha(p/'initial_flow.csv'),sha(ROOT/'outputs/v126_validation/format01/restart_flow.csv'))
                continue
            with (p/previous).open() as f:oldkeys=next(csv.reader(f))
            old=np.loadtxt(str(p/previous),delimiter=',',skiprows=1)
            with (p/'restart_flow.csv').open() as f:keys=next(csv.reader(f))
            a=np.loadtxt(str(p/'restart_flow.csv'),delimiter=',',skiprows=1)
            self.assertEqual(keys,oldkeys);rho=a[:,keys.index('Density')]
            uv=a[:,[keys.index('Momentum_x'),keys.index('Momentum_y')]]/rho[:,None]
            thermal=a[:,keys.index('Energy')]/rho-np.sum(uv*uv,axis=1)/2-old[:,keys.index('Turb_Kin_Energy')]
            c=config(p/'case.cfg');factor=(float(c['GAMMA_VALUE'])-1)/float(c['GAS_CONSTANT'])
            expected=thermal*factor;native=a[:,keys.index('Temperature')]
            error=float(max(abs(expected/native-1)));self.assertLess(error,2e-12);errors.append(error)
            current=expected+(old[:,keys.index('Turb_Kin_Energy')]-a[:,keys.index('Turb_Kin_Energy')])*factor
            changes.append(float(max(abs(native-current))));self.assertGreater(changes[-1],1e-8)
        self.__class__.metrics=dict(previousKTemperatureRelativeErrors=errors,currentKTemperatureDifferenceK=changes)
    def test_expanded_fields_need_explicit_contract(self):
        p=ROOT/'outputs/v126_validation/format01/restart_flow.csv'
        with self.assertRaises(AssertionError):state(p,True)
        self.assertEqual(state(p,True,True)[2],2)
    def test_missing_or_reordered_field_rejected(self):
        p=ROOT/'outputs/v126_validation/format01/restart_flow.csv';lines=p.read_text().splitlines()[:3]
        with tempfile.TemporaryDirectory() as directory:
            q=Path(directory)/'state.csv'
            for name in ['Temperature','Omega']:
                q.write_text('\n'.join([lines[0].replace(name,'UnknownField')]+lines[1:]))
                with self.assertRaises(AssertionError):state(q,True,True)
    def test_changed_physics_source_rejected(self):
        p=ROOT/'outputs/v126_validation/restart02';c=config(p/'case.cfg');init=read(p/'sst_initialization.json')
        for key,value in [('PRANDTL_TURB','.8'),('CFL_NUMBER','.25'),('MACH_NUMBER','5')]:
            with self.subTest(key=key),self.assertRaises(AssertionError):validate_source(p,dict(c,**{key:value}),init)

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(OutputTests);result=unittest.TextTestRunner().run(suite)
    p=ROOT/'outputs/v126_validation/output_contract_tests.json';assert not p.exists()
    p.write_text(json.dumps(dict(passed=result.wasSuccessful(),tests=result.testsRun,metrics=getattr(OutputTests,'metrics',{}),
        sourceSHA256=sha(Path(__file__).resolve())),indent=2)+'\n',encoding='utf-8')
    raise SystemExit(0 if result.wasSuccessful() else 1)
