"""模块用途。将真实异常地址映射到固定二进制符号，并独立核对同源首帧和退出身份。"""
import bisect
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from auditSU2IdealAirState import load_state
from su2ThermalRuntime import read, sha
from validateIdealAir import save


def main():
    folder = ROOT/'outputs/v163_validation/real02'
    output = ROOT/'outputs/v163_validation/native_review.json'
    assert not output.exists()
    receipt = read(folder/'engine_process.json')
    runtime = read(ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-observed-v1/runtime.json')
    executable = Path(runtime['executable'])
    assert sha(executable) == receipt['inputSHA256'][str(executable)]
    assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
    log = (folder/'engine_stdout.log').read_text('utf8')
    assert log.count('RAH_NATIVE_TERMINATION_BEGIN') == log.count('RAH_NATIVE_TERMINATION_END') == 1
    assert 'what=inverse target outside temperature range' in log
    tool = ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin'
    header = subprocess.check_output([str(tool/'objdump.exe'), '-p', str(executable)], timeout=30).decode()
    image_base = int(re.search(r'ImageBase\s+([0-9a-fA-F]+)', header)[1], 16)
    command = [str(tool/'nm.exe'), '-n', '-C', '--defined-only', str(executable)]
    text = subprocess.check_output(command, timeout=30).decode('utf8', errors='replace')
    symbols = []
    for line in text.splitlines():
        match = re.fullmatch(r'([0-9a-fA-F]+) [TtWw] (.*)', line)
        if match:
            symbols.append((int(match[1], 16), match[2]))
    symbols.sort()
    addresses = [value[0] for value in symbols]
    frames = []
    for match in re.finditer(r'frame=(\d+) address=([0-9a-fA-F]+) module=(.*?) base=([0-9a-fA-F]+) offset=0x([0-9a-fA-F]+)', log):
        if Path(match[3]) != executable:
            continue
        address = image_base+int(match[5], 16)
        symbol_address, name = symbols[bisect.bisect_right(addresses, address-1)-1]
        frames.append(dict(frame=int(match[1]), moduleRelativeAddress=match[5],
            preferredAddress=hex(address), symbol=name, symbolOffset=address-symbol_address))
    assert frames
    old = ROOT/'outputs/v161_validation/control'
    first, fields, nd = load_state(folder/'restart_flow_000000.csv', True, captured=True)
    previous, old_fields, old_nd = load_state(old/'restart_flow_000000.csv', True, captured=True)
    assert list(fields) == list(old_fields) and first.shape == previous.shape and nd == old_nd == 3
    difference = float(np.max(abs(first-previous)/np.maximum(1, abs(previous))))
    assert difference < 2e-13
    report = dict(exceptionType=re.search(r'exception_type=(.*)', log)[1],
        exceptionWhat=re.search(r'what=(.*)', log)[1], stack=frames,
        nativeExitCode=receipt['exitCode'], nativeSeconds=receipt['elapsedSeconds'],
        inputsStable=all(receipt['inputsStable'].values()), lockedInputs=len(receipt['inputsStable']),
        checkpoints=receipt['checkpoints'], firstFrameRelativeDifference=difference,
        symbolCommand=command, preferredImageBase=hex(image_base),
        physicalCauseQualified=False, physicsQualified=False, newQualifiedArea=0,
        inputSHA256={str(p): sha(p) for p in [executable, folder/'engine_process.json',
            folder/'engine_stdout.log', old/'restart_flow_000000.csv', folder/'restart_flow_000000.csv',
            Path(__file__).resolve(), tool/'nm.exe', tool/'objdump.exe']})
    save(output, report)
    print(json.dumps(dict(exception=report['exceptionWhat'], frames=frames,
        nativeSeconds=report['nativeSeconds'], firstFrameRelativeDifference=difference), indent=2))


if __name__ == '__main__':
    main()
