"""空间诊断的实际数据、严格拒绝和解析中位体积校核；不启动流场。"""
from pathlib import Path
import sys,json,csv,unittest,tempfile,copy
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import sha,read,config
from su2SSTInterface import state
from su2ThermalRuntime import validate_sst_manifest
from diagnoseSU2SSTSpatial import median_volumes
from prepareSU2SSTCheckpoint import validate_source


class SpatialTests(unittest.TestCase):
    folder=ROOT/'outputs/v128_validation/spatial01'

    def test_actual_terminal_and_rms(self):
        r=read(self.folder/'engine_process.json');s=read(self.folder.parent/'spatial_analysis.json')
        self.assertEqual(r['status'],'terminal');self.assertEqual(r['exitCode'],0)
        self.assertTrue(all(r['inputsStable'].values()));self.assertEqual(s['nodes'],349206)
        self.assertLess(s['maximumLogRMSDifference'],2e-8);self.assertEqual(s['originalFaces'],19752)
        self.assertFalse(s['physicsQualified']);self.assertFalse(s['originalResidualGatePassed'])
        self.assertLess(abs(sum(g['volumeM3'] for g in s['layerGroups'].values())/sum(g['volumeM3'] for g in s['cellTypes'].values())-1),2e-10)

    def test_new_columns_require_explicit_identity(self):
        with (self.folder/'restart_flow.csv').open() as f:rows=[next(f) for _ in range(3)]
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'state.csv';p.write_text(''.join(rows))
            with self.assertRaises(AssertionError):state(p,True,True)
            self.assertEqual(state(p,True,True,True)[2],3)
            rows[0]=rows[0].replace('Residual_Energy','Undeclared_Energy')
            p.write_text(''.join(rows))
            with self.assertRaises(AssertionError):state(p,True,True,True)

    def test_analytic_prism_and_tetra_median_volumes(self):
        p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,0,1],[0,1,1],
                    [2,0,0],[3,0,0],[2,1,0],[2,0,1]],float)
        v,q,total=median_volumes(p,[[13,0,1,2,3,4,5],[10,6,7,8,9]])
        np.testing.assert_allclose(v[:6],1/12,rtol=2e-15)
        np.testing.assert_allclose(v[6:],1/24,rtol=2e-15)
        self.assertAlmostEqual(v.sum(),2/3)

    def test_same_configuration_protection_retained(self):
        c=config(self.folder/'case.cfg');init=read(self.folder/'sst_initialization.json')
        for key,value in [('CFL_NUMBER','.2'),('CONV_NUM_METHOD_FLOW','ROE'),('PRANDTL_TURB','.8')]:
            with self.subTest(key=key),self.assertRaises(AssertionError):
                validate_source(self.folder,dict(c,**{key:value}),init)

    def test_diagnostic_cannot_be_silently_extended(self):
        # 配置摘要同时更新，仍必须由1步诊断保护拒绝，而不是仅被过期SHA拒绝。
        text=(self.folder/'case.cfg').read_text();data=read(self.folder/'case_manifest.json')
        for old,new in [('ITER= 1','ITER= 2'),('VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL, TIMESTEP)',
                                         'VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL)')]:
            with self.subTest(change=new),tempfile.TemporaryDirectory() as d:
                p=Path(d);cfg=text.replace(old,new);self.assertNotEqual(text,cfg)
                (p/'case.cfg').write_text(cfg);(p/'mesh.su2').write_text('negative-test mesh')
                m=copy.deepcopy(data);m['sst']['configurationSHA256']=sha(p/'case.cfg');m['sst']['meshSHA256']=sha(p/'mesh.su2')
                (p/'sst_contract.json').write_text(json.dumps(m['sst']));m['sstContractSHA256']=sha(p/'sst_contract.json')
                with self.assertRaises(AssertionError):validate_sst_manifest(p,m,cfg)


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(SpatialTests);r=unittest.TextTestRunner().run(suite)
    dest=ROOT/'outputs/v128_validation/spatial_tests.json';assert not dest.exists()
    dest.write_text(json.dumps(dict(passed=r.wasSuccessful(),tests=r.testsRun,sourceSHA256=sha(Path(__file__).resolve()),
        sourceBindings={str(p):sha(p) for p in [ROOT/'scripts/su2SSTInterface.py',ROOT/'scripts/su2ThermalRuntime.py',
            ROOT/'scripts/diagnoseSU2SSTSpatial.py',ROOT/'scripts/prepareSU2SSTCheckpoint.py']}),indent=2)+'\n',encoding='utf-8')
    raise SystemExit(0 if r.wasSuccessful() else 1)
