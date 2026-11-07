"""固定原生热物性分派与Cp积分源码审核；不构建、不运行流场。"""
from pathlib import Path
import json,hashlib,shutil,time
ROOT=Path(__file__).resolve().parents[1]


def main():
    start=time.monotonic();out=ROOT/'outputs/v202_validation/source_review02';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');engine=Path('/opt/openfoam13')
    names=['src/thermophysicalModels/specie/include/forGases.H','src/thermophysicalModels/specie/include/forThermo.H',
        'src/thermophysicalModels/basic/psiThermo/psiThermos.C',
        'src/thermophysicalModels/specie/thermo/hIcoTabulated/hIcoTabulatedThermo.H',
        'src/thermophysicalModels/specie/thermo/hIcoTabulated/hIcoTabulatedThermo.C',
        'src/thermophysicalModels/specie/thermo/hIcoTabulated/hIcoTabulatedThermoI.H',
        'src/thermophysicalModels/specie/thermo/hTabulated/hTabulatedThermoI.H',
        'src/thermophysicalModels/specie/equationOfState/perfectGas/perfectGasI.H']
    for name in ['makeFluidThermo.H','thermodynamicConstants.C']:
        matches=[p.relative_to(engine).as_posix() for p in engine.glob('src/**/'+name) if 'lnInclude' not in p.parts]
        if len(matches)!=1:raise ValueError('必须唯一准确源 '+name)
        names+=matches
    for pattern in ['src/**/integratedNonUniformTable1.C','src/**/integratedNonUniformTable1.H',
        'src/**/nonUniformTable1.C','src/**/nonUniformTable1.H','src/**/nonUniformTable1I.H']:
        names += [p.relative_to(engine).as_posix() for p in engine.glob(pattern) if 'lnInclude' not in p.parts]
    rows=[]
    for n in names:
        p=engine/n
        if not p.is_file():raise ValueError('准确路径缺失 '+n)
        q=out/'source'/n;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=n,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
    text=(engine/names[0]).read_text();ico=(engine/names[5]).read_text();tab=(engine/names[6]).read_text()
    assert 'hIcoTabulated' not in text and 'Cp_.intfdT(T)' in ico and 'Cp_.intfByTdT(T)' in ico and 'NotImplemented' in tab
    r=dict(passed=True,stockGasRegistration=False,selectedNativeClass='hIcoTabulatedThermo<perfectGas<specie>>',
        requiredAdapter='register native constTransport / sensibleInternalEnergy / psiThermo / pureMixture combination',
        rejectedAlternative='hTabulated s and gStd are NotImplemented; independent H/Cp interpolation',
        files=rows,engineBuilds=0,CFDRuns=0,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))


if __name__=='__main__':main()
