"""用九个已完成阶段的回执检查恢复行为。数值场的独立复算由 verify_later 负责。"""
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
import shortSupportLater65 as solver


def check(job):
    work=job.resolve()/"run"
    saved=solver.read(job/"acceptedSteps.json")["steps"]
    assert len(saved)==9
    def inspect(work,branch,k,*args,**kwargs):
        return work/branch/f"step{k}"/f"restart_{k:05d}.csv"
    def energy(work,branch,k,**kwargs):
        return solver.read(work/branch/f"step{k}"/"energy.json")
    with tempfile.TemporaryDirectory() as tmp:
        test=Path(tmp)
        journal=test/"acceptedSteps.json"
        for length,damaged in ((9,False),(8,False),(9,True)):
            rows=[dict(row) for row in saved[:length]]
            if damaged:rows[3]["old"]="different-source"
            journal.write_text(json.dumps(dict(steps=rows)),encoding="utf8")
            original=journal.read_bytes()
            # 使用真实能量账和源帧指纹，只拦住另行验收的原生求解与重复读场。
            with patch.object(solver,"advance",side_effect=AssertionError("已存步骤被重新启动")), \
                 patch.object(solver,"inspect_step",side_effect=inspect),patch.object(solver,"endpoint",side_effect=energy), \
                 patch.object(solver,"progress"),patch.object(solver,"atomic_write_json",wraps=solver.atomic_write_json) as writes:
                if damaged:
                    try:solver.advance_chain(test,work)
                    except ValueError as err:assert "能量账发生变化" in str(err)
                    else:raise AssertionError("错误的旧回执被接受")
                else:
                    assert solver.advance_chain(test,work)==saved
                assert writes.call_count==int(length<len(saved))
                if length==len(saved):assert journal.read_bytes()==original
                else:assert solver.read(journal)["steps"]==saved
    print("Checked intact, partial and damaged ledgers. No accepted physical step was restarted.")


if __name__=="__main__":check(Path(sys.argv[1]))
