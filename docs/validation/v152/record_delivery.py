"""保存本轮小型证据和逐字执行源，完整网格/场继续保留在原运行目录。"""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OUT = ROOT/'outputs/v152_validation'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def copy_exact(source, target):
    data = source.read_bytes()
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        assert target.read_bytes() == data, target
    else:
        target.write_bytes(data)


sources = ['scripts/prepareSU2AxisymmetricContinuation.py', 'scripts/su2AxisymmetricSurface.py',
           'scripts/runSU2AxisymmetricStudy.py', 'scripts/auditSU2AxisymmetricStudy.py',
           'scripts/su2AxisymmetricAnalysis.py', 'scripts/su2BenchmarkAcceptance.py',
           'scripts/prepareSU2AxisymmetricInitialState.py',
           'readSU2BenchmarkReview.m', 'createSU2BenchmarkReview.m',
           'tests/testSU2AxisymmetricContinuation.py', 'tests/testSU2BenchmarkContinuationReview.m']
for name in sources:
    copy_exact(ROOT/name, HERE/'executed_sources'/name)
names = ['small_test02.log', 'final_small_test.log', 'prepare02.log', 'native_run.log', 'audit.log',
         'delivery_verification02.log', 'saved_diagnostics.log', 'ui_process.json', 'ui_stdout.log',
         'review01/report.json', 'review01/convergence.png', 'workbench01/report.json',
         'workbench01/benchmark_review.json.stations.csv',
         'acc_sbli01/study_manifest.json', 'acc_sbli01/engine_process.json',
         'acc_sbli01/initialization_run_binding.json', 'acc_sbli01/continuation.json',
         'acc_sbli01/surface_pairing.json', 'acc_sbli01/case.cfg', 'acc_sbli01/history.csv',
         'acc_sbli01/engine.log', 'acc_sbli01/surface_flow_005700.vtk']
for name in names:
    copy_exact(OUT/name, HERE/'evidence'/name)

# 新路线只记录准确原始理想项及EOS元数据，实际有效温区和实现仍待审核。
air = ROOT/'tools/runtime/su2-thermal-build/source/subprojects/CoolProp/dev/fluids/Air.json'
eos = json.loads(air.read_text())['EOS'][0]
review = dict(source=str(air), sourceSHA256=sha(air),
              sourceCommit='98b3523d5daa98454618d381d2ae53f7471d216b',
              data={key: eos[key] for key in ('BibTeX_CP0', 'BibTeX_EOS', 'T_max', 'Ttriple',
                                             'gas_constant', 'molar_mass', 'alpha0')},
              implemented=False, tested=False, independentIdealTermTemperatureRangeVerified=False)
target = HERE/'source_review'
target.mkdir(exist_ok=True)
(target/'air_ideal_candidate.json').write_text(json.dumps(review, indent=2)+'\n', encoding='utf8')
archive = ROOT/'tools/runtime/su2-thermal-build/su2.zip'
z = zipfile.ZipFile(str(archive))
excerpts = []
for suffix, term, after in [('Common/src/CConfig.cpp', 'addDoubleOption("MAX_UPDATE_FLOW"', 7),
                          ('SU2_CFD/src/solvers/CEulerSolver.cpp', 'void CEulerSolver::ComputeUnderRelaxationFactor', 35)]:
    name = next(n for n in z.namelist() if n.endswith(suffix))
    lines = z.read(name).decode().splitlines()
    index = next(i for i, line in enumerate(lines) if term in line)
    excerpts.append(dict(member=name, firstLine=index+1, lines=lines[index:index+after]))
(target/'native_update_source.json').write_text(json.dumps(dict(archiveSHA256=sha(archive),
    excerpts=excerpts, scope='source review only; NaN cause unproven; no configuration change'), indent=2)+'\n', encoding='utf8')
manifest = dict(executedSources={name: sha(ROOT/name) for name in sources},
                evidence={name: sha(OUT/name) for name in names},
                uniqueNumericalRuns=1, numericalExitCode=1, matlabRuns=1, matlabExitCode=0,
                newQualifiedArea=0, newNativeMouseOrScreenshot=False,
                largeFieldsLocalOnly=True, uniqueDataDeleted=False, projectArchiveCreated=False)
(HERE/'delivery_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(sourceFiles=len(sources), evidenceFiles=len(names))))
