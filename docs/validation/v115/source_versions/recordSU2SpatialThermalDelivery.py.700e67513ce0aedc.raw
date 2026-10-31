"""归档v115真实运行/构建来源；只读科学输出，不重算、不把失败改成成功。"""
from pathlib import Path
import argparse,hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v115_validation';DEST=ROOT/'docs/validation/v115'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def copy(p,q):
    q.parent.mkdir(parents=True,exist_ok=True)
    if q.exists():assert sha(p)==sha(q),str(q)
    else:shutil.copyfile(str(p),str(q))
def source(p):
    h=sha(p);q=DEST/'source_versions'/(p.name+'.'+h[:16]+'.raw');copy(p,q)
    return dict(sha256=h,archived=q.relative_to(ROOT).as_posix())
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--finalize',action='store_true');args=parser.parse_args()
    assert not (DEST/'delivery_manifest.json').exists()
    info=read(OUT/'actual01/report.json');engine=Path(info['engineFolder'])
    cases={'official_constant':OUT/'actual01/official_constant','native_constant':OUT/'actual01/native_constant','spatial':engine}
    for name,p in cases.items():
        receipt=read(p/'engine_process.json');assert receipt['status']=='terminal' and receipt['exitCode']==0
        assert all(receipt['inputsStable'].values())
        for n in ['case.cfg','case_manifest.json','engine_request.json','engine_process.json','engine_stdout.log','supervisor.log',
                  'history.csv','surface_flow.vtk','wall_temperature.csv','rah_wall_temperature_applied.csv']:
            if (p/n).exists():copy(p/n,DEST/'cases'/name/n)
    for p in (OUT/'native_rejections').glob('*/*'):
        if p.name in ['process.json','engine_stdout.log','case.cfg','wall_temperature.csv']:
            copy(p,DEST/'rejections'/p.parent.name/p.name)
    for p in OUT.iterdir():
        if p.is_file() and p.suffix in ['.json','.log','.diff']:copy(p,DEST/'processes'/p.name)
    for part in ['actual01','review01','unit01']:
        for p in (OUT/part).glob('*'):
            if p.is_file() and p.suffix in ['.json','.csv','.mat']:
                copy(p,DEST/'checks'/part/p.name)
    build=ROOT/'tools/runtime/su2-thermal-build'
    for p in build.glob('*_download.json'):copy(p,DEST/'build'/p.name)
    for rel in ['build/meson-info/intro-compilers.json','build/meson-info/intro-buildoptions.json',
                'build/meson-logs/meson-log.txt','build/build.ninja','build/compile_commands.json']:
        p=build/rel
        if p.exists():copy(p,DEST/'build'/p.name)
    # 准确封存本轮构建触及的上游版本；官方许可证沿用third_party/su2/COPYING。
    compiled=read(OUT/'build_sources04.json')
    for path,digest in compiled.items():
        p=Path(path)
        if p.suffix in ['.cpp','.hpp','.inc','.py']:
            assert sha(p)==digest;source(p)
    used=['aeroHeatWorkflowState.m','buildSchaeferFergusonBody.m','schaeferFergusonCondition.m',
          'resolveSU2ThermalBoundary.m','attachSU2ThermalBoundary.m','writeSU2RocketCase.m',
          'readSU2HybridVolume.m','su2WorkbenchMeshOptions.m','rocketAeroHeatWorkbench.m',
          'createSU2ResearchView.m','createAerothermalGeometryView.m','updateAerothermalGeometryView.m',
          'buildSU2WallProjection.m','readSU2SurfaceResult.m','mapSU2WallResult.m','embeddedSU2Job.m',
          'auditSU2ResultConvergence.m','describeSU2Convergence.m','tests/testSU2SpatialThermalWorkbench.m',
          'tests/reviewSU2SpatialBoundaryUI.m','tests/testSchaeferThermalContract.m']
    used += ['web/'+n for n in ['geometry-su2.html','geometry.html','geometry.css','geometry.js','ale-section-core.js','su2-section-view.js']]
    jobs=[]
    for name in ['spatial01','review01']:
        p=OUT/(name+'_process.json');r=read(p)
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['sourceStableAfter'].values())
        exact={}
        for rel in used:
            if rel not in r['sourceSHA256']:continue
            assert sha(ROOT/rel)==r['sourceSHA256'][rel];exact[rel]=source(ROOT/rel)
        jobs.append(dict(name=name,exitCode=0,sourceStable=True,sources=exact))
    if not args.finalize:
        print('Exact build and terminal numerical/UI evidence archived before cache cleanup.');return
    sources={rel:source(ROOT/rel) for rel in used}
    for rel in ['scripts/runEmbeddedSU2.py','scripts/su2ThermalRuntime.py','scripts/buildSU2ThermalRuntime.py',
                'scripts/prepareSU2ThermalBuild.py','scripts/installSU2ThermalRuntime.py','scripts/recordSU2SpatialThermalDelivery.py',
                'tests/auditSU2V115.py','tests/auditSU2V110.py','tests/testSU2ThermalRuntimeRejection.py','tests/testSU2AxialThermalParser.cpp',
                'third_party/su2_thermal/RAHAxialWallTemperature.hpp','third_party/su2_thermal/SetCustomBoundary.cpp.inc']:
        sources[rel]=source(ROOT/rel)
    manifest=dict(baseCommit='7b05564c7f761b5d1f4d0ec4a969da14429a3c0d',
        status='tested-spatial-thermal-interface; three-dimensional-flow-unconverged',matlabJobs=jobs,
        sources=sources,newQualifiedArea=0,physicsQualified=False,
        evidence={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(DEST.rglob('*')) if p.is_file()},
        localOnly=['fixed official and controlled runtime binaries','full mesh and restart fields','shared historical store'])
    (DEST/'delivery_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(manifest['evidence']),sources=len(sources),matlabTerminalJobs=len(jobs))))
if __name__=='__main__':main()
