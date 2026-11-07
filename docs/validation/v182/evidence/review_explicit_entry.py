"""定向保存固定引擎的显式入口和步长限制；只读源码，不构建或运行CFD。"""
from pathlib import Path
import hashlib
import json
import subprocess
import time

start = time.monotonic()
root = Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source')
out = Path(__file__).resolve().parent / 'explicit_entry01'
out.mkdir(exist_ok=False)
revision = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
assert revision == 'f53f4609a0331d48efee69a4e4f3c3598378cc03'
assert not subprocess.check_output(['git', '-C', str(root), 'diff', '--name-only'], text=True).strip()

# 模块一：只取已定位的准确区段，并绑定完整文件SHA，不以文件名猜测行为。
spans = {
    'globalconfig.d': [(103, 116), (154, 173), (220, 235), (2418, 2430)],
    'lua-modules/globalconfig.lua': [(10, 75), (125, 160)],
    'fluidfvcell.d': [(998, 1070)],
    'fluidblock.d': [(1182, 1218)],
    'simcore_gasdynamic_step.d': [(1088, 1208)],
    'timemarching.d': [(420, 448)],
    'commands/runsim.d': [(342, 368)],
}
files = []
for rel, ranges in spans.items():
    src = root / 'src/lmr' / rel
    raw = src.read_bytes()
    lines = raw.decode('utf8').splitlines()
    excerpt = '\n\n'.join('\n'.join('%d: %s' % (i, lines[i - 1]) for i in range(lo, hi + 1)) for lo, hi in ranges) + '\n'
    dest = out / (rel.replace('/', '__') + '.txt')
    dest.write_text(excerpt, encoding='utf8')
    files.append(dict(source=str(src), sha256=hashlib.sha256(raw).hexdigest(), ranges=ranges,
                      excerpt=dest.name, excerptSHA256=hashlib.sha256(dest.read_bytes()).hexdigest()))

# 模块二：这些是准确源码职责，不是实际显式运行、内存或物理验证回执。
report = dict(sourceReviewCompleted=True, sourceCommit=revision, files=files,
              proposedProfile='rah-explicit-local-v1', implemented=False, nativeRuns=0, newCFDSteps=0,
              confirmed=['time-marching dispatches to native time integration',
                         'tvd_rk3 has three update stages and four flow time levels',
                         'local dt uses acoustic/convective and viscous signals plus turbulent source restriction',
                         'dt capped by global dt times local factor and dt_max',
                         'integrateInTime failure may only print a message; exit0 is insufficient'],
              requiredInputs=['fixed native lmr-run and original thermally-perfect air/Wilcox model',
                              'original mesh/wall mapping, Tw, freestream, turbulence and boundary data',
                              'explicit RK3/local stepping/stringent CFL and fixed diffusion/turbulence factors',
                              'interpolation/flux/limiter, initial local dt and local ratio cap',
                              'max step/pseudo-time/wall-clock, field/load cadence, native completion audit',
                              'group memory/disk/CPU limits and cleanup receipts'],
              limitations=['pseudo-time is not a physical trajectory',
                           'no proof of full-domain memory fit or convergence',
                           'thin near-wall cells may impose severe diffusion stiffness',
                           'old NK iteration comparison is not an explicit integration acceptance test'],
              physicalQualification=False, newQualifiedAreaM2=0, elapsedSeconds=time.monotonic()-start)
(out / 'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
print(json.dumps(dict(sourceReviewCompleted=True, files=len(files), elapsedSeconds=report['elapsedSeconds'])))
