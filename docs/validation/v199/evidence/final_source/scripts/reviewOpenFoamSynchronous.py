"""绑定成熟同步Euler/upwind组合的准确源；只读选择，不重编或求解。"""
from pathlib import Path
import json,time,shutil
from runEilmerProcess import sha

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v199_validation'
FOAM=Path('/opt/openfoam13')

def main():
    start=time.monotonic();out=BASE/'source_review01';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    pkg=json.loads((ROOT/'tools/runtime/openfoam/13-20260624/runtime.json').read_text())
    for path,h in pkg['filesSHA256'].items():assert sha(path)==h
    extra=['applications/solvers/foamRun/foamRun.C',
        'src/finiteVolume/finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.C',
        'src/finiteVolume/finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.C',
        'src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H',
        'src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/vanLeer/vanLeer.H',
        'src/OpenFOAM/matrices/lduMatrix/solvers/diagonalSolver/diagonalSolver.C',
        'tutorials/shockFluid/wedge15Ma5/system/fvSchemes',
        'tutorials/shockFluid/wedge15Ma5/system/controlDict']
    sources=[FOAM/p for p in extra]+list((FOAM/'applications/modules/shockFluid').glob('*.C'))
    for p in sources:
        q=out/'source'/p.relative_to(FOAM);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    fluid=(FOAM/'applications/solvers/foamRun/foamRun.C').read_text()
    assert fluid.index('solver.prePredictor();')<fluid.index('solver.momentumPredictor();')<fluid.index('solver.thermophysicalPredictor();')<fluid.index('solver.pressureCorrector();')
    upwind=(FOAM/extra[3]).read_text();assert 'return pos0(this->faceFlux_);' in upwind
    euler=(FOAM/extra[1]).read_text();assert 'const scalar rDeltaT = 1.0/mesh().time().deltaTValue();' in euler
    energy=(FOAM/'applications/modules/shockFluid/thermophysicalPredictor.C').read_text()
    assert 'fvm::ddt(rho, e)' in energy and 'fvc::ddt(rho, K)' in energy and 'thermo_.correct();' in energy
    report=dict(passed=True,profile='of13-sync-euler-upwind-nasa7-v1',
        actualChanges=['one global dt replaces per-cell pseudo-time','piecewise-constant complete face states replace separate limited reconstructions'],
        numericalOrder='first-order; accuracy/refinement remain required',generalVariableCpPositivityTheorem=False,
        sourceSHA256={str(p):sha(p) for p in sources},runtimeSHA256=sha(ROOT/'tools/runtime/openfoam/13-20260624/runtime.json'),
        paperURL='https://www.math.umd.edu/~tadmor/centpack/publications/files/Greenshields%20et%20al%20polyhedral%20FV%20IJNMF2009-centpack.pdf',
        paperScope='calorically perfect only; no NASA7 qualification',
        CFDRuns=0,engineBuilds=0,elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(passed=True,elapsedSeconds=report['elapsedSeconds'])))

if __name__=='__main__':main()
