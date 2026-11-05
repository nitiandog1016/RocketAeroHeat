"""复用v121几何执行两条实际非零热流路径，有限更新仅取得接口证据。"""
from pathlib import Path
import os
import shutil
import sys

from validateIdealAir import ROOT, run, save, sha
from prepareSU2SSTCase import set_value
from prepareSU2IdealAirCase import prepare


def main():
    source = ROOT/'outputs/v121_validation/capability01'
    out = ROOT/'outputs/v154_validation/cases01'
    out.mkdir(exist_ok=False)
    # 两组各10更新、逐次保存相邻k阶段，预计场文件总长度约50MiB。
    # 每组60秒、两线程；外层监督另外限制总工作集4GiB。
    for name, turbulent in (('laminar', False), ('sst', True)):
        folder = out/name
        folder.mkdir()
        os.link(str(source/'mesh.su2'), str(folder/'mesh.su2'))
        text = (source/'case.cfg').read_text(encoding='utf8')
        text = set_value(text, 'ITER', '10')
        (folder/'case.cfg').write_text(text, encoding='utf8')
        save(folder/'case_manifest.json', dict(schemaVersion=1, sourceGeometry=str(source/'mesh.su2'),
            sourceMeshSHA256=sha(source/'mesh.su2'), physicsQualified=False))
        prepare(folder, turbulent, 'v154 finite nonzero-wall-heat interface check; assumed turbulence parameters')
        sources = ['scripts/prepareSU2IdealAirCase.py', 'scripts/su2IdealAirInterface.py',
                   'scripts/su2SSTInterface.py', 'scripts/validateSU2IdealAirCases.py',
                   'docs/validation/v154/PLAN.md']
        save(folder/'engine_request.json', dict(workDirectory=str(folder), configFile='case.cfg',
            meshFile='mesh.su2', maximumSeconds=60, threads=2, dryRun=False, sourceFiles=sources))
        snapshots = out/(name+'_executed_sources')
        snapshots.mkdir()
        for filename in sources+['scripts/runEmbeddedSU2.py', 'scripts/su2ThermalRuntime.py']:
            target = snapshots/filename
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(str(ROOT/filename), str(target))
        run(out, name, [sys.executable, '-X', 'utf8', str(ROOT/'scripts/runEmbeddedSU2.py'),
                        '--request', str(folder/'engine_request.json')], 70)
    total = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    save(out/'storage.json', dict(fileLengthBytes=total, limitBytes=80*1024**2, physicalAllocatedBytes=None))
    assert total <= 80*1024**2


if __name__ == '__main__':
    main()
