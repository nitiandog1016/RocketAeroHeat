"""在独立安装中核对缓存，并用临时副本检查源码变化和程序损坏。"""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch


def check(app):
    sys.path[:0] = [str(app/"scripts"), str(app/"docs/validation/v333")]
    import shortSupportTime65
    import shortSupportFirstRun as first
    import buildShort65Programs as programs
    import installShort65Kernel as installer
    from shortSphere65 import read, sha
    assert shortSupportTime65.ROOT == first.ROOT == programs.ROOT == app
    kernel_file=installer.HOME_DIR/"kernel.json"
    before=sha(kernel_file)
    assert installer.install()==read(kernel_file) and sha(kernel_file)==before
    with tempfile.TemporaryDirectory(dir=app/"outputs") as folder:
        scratch = Path(folder)
        preps = {n:sha(first.CACHE/(n+".json")) for n in ("inputs","programs","base")}
        prep, native, steady, stages = first.inputs(scratch)
        assert preps == {n:sha(first.CACHE/(n+".json")) for n in preps}
        assert not (scratch/"progress.json").exists()
        # 单独验证带预算的正式首次建网入口。临时网在确认六个输入相同后释放。
        grid=scratch/"prepared"
        first.prepare_grid(grid)
        mesh_run=read(scratch/"mesh_process.json")
        assert not (mesh_run["stopReason"] or mesh_run["livePids"] or mesh_run["exitCode"])
        assert sha(grid/"prepared.json")==sha(prep/"prepared.json")
        assert all(sha(grid/name)==key for name,key in read(prep/"prepared.json")["files"].items())
        # 替换生成器中的源片段，核对生产程序拒绝与本版本不相符的缓存。
        with patch.object(programs, "APPLY", programs.APPLY+"\n// changed test fixture\n"):
            try:
                programs.check_programs(native, prep)
            except ValueError as err:
                assert "源码已经改变" in str(err)
                source = str(err)
            else:
                raise AssertionError("变更过的适配源码命中了旧程序")
        # 使用独立空文件模拟损坏，其余文件只读引用。原始可执行文件始终保持完整。
        bad = scratch/"damaged"
        for kind in ("steady","time"):
            (bad/kind).mkdir(parents=True)
            for f in (native/kind).iterdir():
                target = bad/kind/f.name
                if f.name=="physical.exe":
                    target.write_bytes(b"invalid build")
                elif f.is_file():
                    os.link(f,target)
        try:
            programs.check_programs(bad,prep)
        except ValueError as err:
            assert "程序缓存与编译记录不符" in str(err)
            binary = str(err)
        else:
            raise AssertionError("已损坏的程序命中了缓存")
        return dict(cacheBeforeAndAfter=preps, kernelBeforeAndAfter=before, gridProcess=mesh_run,
                    prepared=sha(prep/"prepared.json"), steady=stages["frame"],
                    steadyStage=stages["stage"], thermal=sha(steady/"thermal.json"), sourceRejection=source,
                    executableRejection=binary, finalExecutable=sha(native/"time/physical.exe"),
                    recordedExecutable=read(native/"time/build.json")["files"]["physical.exe"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("app",type=Path);p.add_argument("output",type=Path)
    args=p.parse_args();report=check(args.app.resolve())
    assert report["finalExecutable"]==report["recordedExecutable"]
    with args.output.open("x",encoding="utf8") as stream:
        json.dump(report,stream,ensure_ascii=False,indent=2)
        stream.write("\n")
