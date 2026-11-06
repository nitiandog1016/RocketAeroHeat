"""只读复核真实返回地址属于SetPrimVar首次物性查询，保存固定二进制和准确源证据。"""
from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from su2ThermalRuntime import read, sha
from validateIdealAir import save

base = ROOT/'outputs/v163_validation'
output = base/'callsite_review.json'
assert not output.exists()
review = read(base/'native_review.json')
frame = next(row for row in review['stack'] if row['symbol'].startswith('CNSVariable::SetPrimVar('))
assert frame['preferredAddress'] == '0x1404804a3'
exe = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-observed-v1/bin/SU2_CFD.exe'
tool = ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/objdump.exe'
command = [str(tool), '-d', '-C', '--start-address=0x14047ffa0', '--stop-address=0x140480650', str(exe)]
assembly = subprocess.check_output(command, timeout=30)
(base/'SetPrimVar_disassembly.txt').write_bytes(assembly)
text = assembly.decode()
assert '140480097:' in text and '140480478' in text and '1404804a1:' in text
assert '140480180:' in text and '140480519:' in text
source = base/'build03/executed_sources/SU2_CFD/src/variables/CNSVariable.cpp'
tree = read(base/'build03/source_tree.json')
original = ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/variables/CNSVariable.cpp'
assert sha(source) == tree[str(original)]
code = source.read_text('utf8')
start = code.index('bool CNSVariable::SetPrimVar(')
end = code.index('\nvoid CNSVariable::SetSecondaryVar(', start)
function = code[start:end]
assert function.count('FluidModel->SetTDState_rhoe(density, staticEnergy);') == 2
assert function.index('FluidModel->SetTDState_rhoe') < function.index('bool check_dens') < function.index('Solution(iPoint,iVar) = Solution_Old')
assert 'catch' not in function
(base/'SetPrimVar_exact.cpp').write_text(function, encoding='utf8')
save(output, dict(passed=True, observedReturnAddress=frame['preferredAddress'],
    firstQueryCallAddress='0x1404804a1', firstQueryDispatchAddress='0x140480097',
    oldSolutionCopyAddress='0x140480180', fallbackQueryCallAddress='0x140480519',
    observedPath='first current-state rhoe query before native primitive checks and old-solution fallback',
    retainedUnknowns=['node identity', 'density/staticEnergy/current k at exception', 'lower or upper temperature bound',
                      'old state validity with current turbulent kinetic energy', 'cause of trial-state excursion'],
    selectedNextChange='IDEAL_AIR_VARIABLE_CP domain_error enters native old-solution admissibility fallback; invalid old state terminates with explicit state context',
    inputSHA256={str(p): sha(p) for p in (exe, tool, source, base/'native_review.json', Path(__file__).resolve())},
    engineBuilds=0, nativeRuns=0, physicsQualified=False, newQualifiedArea=0))
print('First current-state property query confirmed; original fallback was not reached.')
