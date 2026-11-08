"""一次提取固定OpenFOAM区域黏性/热壁/边界源码，记录实际可复用项和限制。"""
from pathlib import Path
import json,shutil,time
from prepareOpenFoamControl import ROOT,save,sha

def main():
    start=time.monotonic();out=ROOT/'outputs/v207_validation/source01';out.mkdir(parents=True,exist_ok=False)
    engine=Path('/opt/openfoam13');paths=[]
    groups=['applications/modules/shockFluid','src/thermophysicalModels/specie/transport/icoTabulated',
        'src/ThermophysicalTransportModels/fluid/laminar/Fourier',
        'src/ThermophysicalTransportModels/fluid/turbulence/eddyDiffusivity',
        'src/finiteVolume/fields/fvPatchFields/derived/waveTransmissive',
        'src/finiteVolume/fields/fvPatchFields/derived/freestream',
        'src/finiteVolume/fields/fvPatchFields/derived/inletOutlet',
        'src/finiteVolume/fields/fvPatchFields/derived/supersonicFreestream',
        'src/functionObjects/field/wallHeatFlux',
        'src/MomentumTransportModels/momentumTransportModels/RAS/SpalartAllmaras']
    for group in groups:
        paths.extend(p for p in (engine/group).iterdir() if p.is_file() and p.suffix in ('.C','.H'))
    tutorial=engine/'tutorials/shockFluid/biconic25-55Run35'
    paths.extend(p for p in tutorial.rglob('*') if p.is_file() and p.stat().st_size<20000)
    records=[]
    for p in paths:
        if time.monotonic()-start>60:raise TimeoutError('来源提取预算耗尽')
        rel=p.relative_to(engine);q=out/'source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        records.append(dict(path=str(p),relativePath=rel.as_posix(),sha256=sha(p),bytes=p.stat().st_size))
    def source(n):return (engine/n).read_text()
    checks=dict(viscousMomentum='momentumTransport->divDevTau(U)' in source(groups[0]+'/momentumPredictor.C'),
        viscousEnergy='thermophysicalTransport->divq(e) + fvc::div(devTauDotU)' in source(groups[0]+'/thermophysicalPredictor.C'),
        tableMu='mu_.value(T)' in source(groups[1]+'/icoTabulatedTransportI.H'),
        tableKappa='kappa_.value(T)' in source(groups[1]+'/icoTabulatedTransportI.H'),
        FourierTemperatureFlux='-fvc::interpolate(this->alpha()*thermo.kappa())*fvc::snGrad(thermo.T())' in source(groups[2]+'/Fourier.C'),
        nativeWallHeatNegatesOutwardFlux='wallHeatFluxBf[patchi] = -qBf[patchi]' in source(groups[8]+'/wallHeatFlux.C'),
        waveTransmissiveUsesFixedGamma='gamma_(dict.lookup<scalar>' in source(groups[4]+'/waveTransmissiveFvPatchField.C'),
        freestreamIsInletOutlet='inletOutletFvPatchField<Type>' in source(groups[5]+'/freestreamFvPatchField.C'),
        officialTutorialLaminar='simulationType  laminar' in source('tutorials/shockFluid/biconic25-55Run35/constant/momentumTransport'))
    if not all(checks.values()):raise ValueError('准确源与采用依据不符')
    (out/'source_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    shutil.copyfile(__file__,out/'executed_entry.py')
    size=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    if size>16*1024**2:raise ValueError('来源记录超预算')
    result=dict(passed=True,sourceChecks=checks,sources=records,bytes=size,elapsedSeconds=time.monotonic()-start,
        fixedEngineCommit='18870c24d21c6b982e2cdec27b2f59738cca5f90',license='GPL-3.0-or-later',
        selected='stock shockFluid + icoTabulated viscosity/conductivity + v206 bounded caloric adapter',
        laminarControlNotFullRocketQualification=True,SAIsNotSAEdwards=True,
        nativeWallHeatSignConversion=-1,fullVariableCpCharacteristicBoundaryVerified=False,
        officialTutorialStartsFromDevelopedState=True,officialTutorialIsNotNewValidation=True,
        CFDRuns=0,builds=0,physicalQualification=False)
    save(out/'report.json',result);print(json.dumps(dict(passed=True,sources=len(records),bytes=size,seconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
