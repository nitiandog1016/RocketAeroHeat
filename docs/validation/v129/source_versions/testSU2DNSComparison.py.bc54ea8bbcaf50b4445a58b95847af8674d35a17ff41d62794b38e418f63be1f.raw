"""公开源身份、热流符号与独立测站选择的定向检查。"""
from pathlib import Path
import sys,unittest,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from compareSU2DNSPlate import reference,station_bracket,REFERENCE


class Contracts(unittest.TestCase):
    def test_reference_heat_and_rounding(self):
        r,a=reference()
        self.assertEqual(a.shape,(330,28));self.assertFalse(r['exactInletAvailable'])
        self.assertGreater(r['heatFluxIntoWallWm2'],10000)
        self.assertLess(r['heatFluxIntoWallWm2'],13000)
        self.assertLess(abs(r['MachFromSpeed']/5.84-1),.003)
        self.assertTrue(all(x['heatFluxIntoWallWm2']>0 for x in r['derivativeCrossChecks']))
        self.assertAlmostEqual(r['heatReferenceRoundingFraction'],1/28)

    def test_tampered_reference_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'reference.dat';p.write_bytes(REFERENCE.read_bytes()+b'\n')
            with self.assertRaises(AssertionError):reference(p)

    def test_station_linear_exact(self):
        i,f=station_bracket([.05,.1,.2,.3],[500,1500,2500,3500],2121,.4)
        self.assertEqual(i,1);self.assertAlmostEqual(f,.621)

    def test_no_extrapolation_or_multiple_crossings(self):
        with self.assertRaises(AssertionError):station_bracket([.1,.2,.3],[10,20,30],2121,.4)
        with self.assertRaises(AssertionError):station_bracket([.05,.1,.15,.2,.25],[1,3,1,3,4],2,.4)
        with self.assertRaises(AssertionError):station_bracket([.01,.02,.03],[1,2,3],1.5,1.)


if __name__=='__main__':unittest.main()
