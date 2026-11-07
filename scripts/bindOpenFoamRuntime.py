"""固定发行包、核心源码和实际运行库绑定；不编译或启动CFD。"""
from pathlib import Path
import hashlib,json,subprocess,time

ROOT=Path('/mnt/e/RocketAeroHeat'); FOAM=Path('/opt/openfoam13')
BASE=ROOT/'outputs/v198_validation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    """模块一：比较准确上游blob，再保存实际包/库与关键算法源码清单。"""
    start=time.monotonic(); out=BASE/'runtime_binding01';out.mkdir(exist_ok=False)
    ref=json.loads((BASE/'openfoam_source01/source_binding.json').read_text())
    matches=[]
    for item in ref['files']:
        raw=(FOAM/item['path']).read_bytes()
        blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        assert blob==item['gitBlobSHA1'],item['path']
        matches.append(dict(path=item['path'],blobSHA1=blob,sha256=sha(FOAM/item['path'])))
    platform=FOAM/'platforms/linux64GccDPInt32Opt'
    bins=[platform/'bin'/n for n in ['foamRun','checkMesh','foamPostProcess']]
    libs=[platform/'lib'/n for n in ['libshockFluid.so','libfluidThermophysicalModels.so',
        'libspecie.so','libfiniteVolume.so','libOpenFOAM.so']]
    sources=list((FOAM/'applications/modules/shockFluid').glob('*.C'))+list((FOAM/'applications/modules/shockFluid').glob('*.H'))
    for rel in ['src/thermophysicalModels/specie/include/forGases.H','src/thermophysicalModels/basic/psiThermo/psiThermos.C',
        'src/thermophysicalModels/specie/specie/specieI.H','src/thermophysicalModels/specie/specie/specie.C',
        'src/thermophysicalModels/specie/thermo/janaf/janafThermoI.H','src/thermophysicalModels/specie/thermo/janaf/janafThermo.C',
        'src/thermophysicalModels/specie/thermo/thermo/thermoI.H',
        'src/finiteVolume/fields/fvPatchFields/derived/freestream/freestreamFvPatchField.C',
        'src/finiteVolume/fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchField.C']:
        sources.append(FOAM/rel)
    for p in sources:
        dst=out/'source'/p.relative_to(FOAM);dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(p.read_bytes())
    version=subprocess.check_output(['dpkg-query','-W','-f=${Version}','openfoam13'],text=True)
    assert version=='20260624'
    helptext=subprocess.check_output(['bash','--noprofile','--norc','-c','source /opt/openfoam13/etc/bashrc; exec foamRun -help'],stderr=subprocess.STDOUT,text=True,timeout=30)
    (out/'help.log').write_text(helptext); assert '13-18870c24d21c' in helptext
    receipt=dict(schema='RAH.OpenFOAMRuntime.v1',profile='openfoam13-20260624-stock',version=version,
        sourceCommit=ref['commit'],upstreamExactFiles=matches,installRoot=str(FOAM),
        packageSHA256='75cb48bb7aee5907fb449ff4b520cbbef5e095c867d8c7320a6ca976aca62791',
        filesSHA256={str(p):sha(p) for p in bins+libs+sources},
        license='GPL-3.0-or-later; unmodified official distribution',fullTreeCommitComparison=False,
        sourceComparisonScope='3 exact upstream blobs plus packaged critical source/binary hashes',
        passed=True,elapsedSeconds=time.monotonic()-start,CFDRuns=0,engineBuilds=0)
    (out/'report.json').write_text(json.dumps(receipt,indent=2)+'\n')
    dst=ROOT/'tools/runtime/openfoam/13-20260624';dst.mkdir(parents=True,exist_ok=False)
    (dst/'runtime.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['passed','elapsedSeconds','version']}))

if __name__=='__main__':main()
