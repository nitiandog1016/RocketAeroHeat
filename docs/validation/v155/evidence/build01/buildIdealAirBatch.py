"""编译独立小型物性批处理库，复用固定头文件，记录真实编译与源身份。"""
import argparse
from pathlib import Path

from validateIdealAir import ROOT, BIN, run, save, sha

RUNTIME = ROOT/'tools/runtime/ideal_air/batch-v1'


def build(output):
    """模块一。只创建新运行库目录，已有构建与失败记录均保留。"""
    output = Path(output).resolve()
    output.relative_to(ROOT/'outputs')
    output.mkdir(parents=True, exist_ok=False)
    RUNTIME.mkdir(parents=True, exist_ok=False)
    sources = [ROOT/'third_party/ideal_air/IdealAir.hpp',
               ROOT/'third_party/ideal_air/IdealAirBatch.cpp', Path(__file__).resolve()]
    locks = {str(path): sha(path) for path in sources}
    save(output/'source_lock.json', locks)
    for path in sources:
        (output/path.name).write_bytes(path.read_bytes())
    library = RUNTIME/'idealAirBatch.dll'
    command = [str(BIN/'g++.exe'), '-std=c++17', '-O2', '-shared',
               '-static-libgcc', '-static-libstdc++', str(sources[1]), '-o', str(library)]
    run(output, 'compile', command, 60)
    assert all(sha(Path(path)) == expected for path, expected in locks.items())
    manifest = dict(schemaVersion=1, model='Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1',
        abi='PT-Ph-rhoe-7fields-v1', library=library.name, sha256=sha(library),
        sourceSHA256=locks, compilerSHA256=sha(BIN/'g++.exe'),
        processSHA256=sha(output/'compile_process.json'))
    save(RUNTIME/'runtime.json', manifest)
    save(output/'runtime.json', manifest)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    build(parser.parse_args().output)
