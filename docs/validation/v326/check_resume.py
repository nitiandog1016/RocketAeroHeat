"""在两个真实拒绝目录检查接续分派，监测场文件、回执及求解调用。"""
from pathlib import Path
import sys
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from shortSupportCase import run, sha


if __name__ == "__main__":
    for version, other in (("v325", "BOUNDED_SCALAR"), ("v326", "SCALAR_UPWIND")):
        job = ROOT/"outputs"/(version+"_validation")/"run01"
        protected = [job/"qualified.json", job/"failure.json", job/"steady01/debugger_process.json",
                     job/"steady01/case.cfg", job/"steady01/history.csv", job/"steady01/restart.csv"]
        before = [sha(p) for p in protected]
        names = sorted(p.name for p in job.iterdir())
        with patch("shortSupportCase.advance", side_effect=AssertionError("不应启动求解")) as advance:
            for method, exception, message in ((other, ValueError, "更换标量对流方式"),
                                               (None, RuntimeError, "没有接受回执")):
                try:
                    run(job, None, None, resume=True, species_convection=method)
                except exception as error:
                    assert message in str(error), error
                else:
                    raise AssertionError("失败阶段被错误接受")
            advance.assert_not_called()
        assert before == [sha(p) for p in protected] and names == sorted(p.name for p in job.iterdir())
        print(version, "resume OK")
