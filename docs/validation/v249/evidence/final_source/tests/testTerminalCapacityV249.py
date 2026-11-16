"""模块职责：快速子进程在首次周期采样前写完超限数据，终态检查仍必须拒绝。"""
from pathlib import Path
import sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import runSU2RocketPhaseV222 as guard
from terminalCapacityV249 import supervise,terminal_capacity


class Capacity(unittest.TestCase):
    def test_fast_successful_process_overflow_is_rejected(self):
        with tempfile.TemporaryDirectory(dir=str(ROOT/'outputs/v249_validation')) as tmp:
            base=Path(tmp);guard.BASE=base;guard.CASE=base;guard.DOC=base
            ok=supervise(guard,base,'fast',[sys.executable,'-c',"from pathlib import Path;Path('burst.dat').write_bytes(b'x'*(2*1024**2))"],10,1,[],base,'offline')
            self.assertFalse(ok);receipt=guard.read(base/'process_fast/report.json');check=guard.read(base/'process_fast/terminal_capacity.json')
            self.assertEqual(receipt['exitCode'],0);self.assertFalse(check['passed']);self.assertGreaterEqual(check['logicalBytes'],2*1024**2)
    def test_outside_process_rejected(self):
        with self.assertRaises(ValueError):terminal_capacity(ROOT/'outputs/v249_validation',ROOT/'outputs/v248_validation/process_rocket',1)


if __name__=='__main__':unittest.main()
