"""模块用途：具名单壁面空间温度层流任务，复用已过网格而不放宽旧小例合同。"""
import argparse
import os
from pathlib import Path
import numpy as np

from naca3344Model import DATA,source_data,wall_temperature,meridian
from prepareSU2NACALaminar import small_geometry,write_mesh,initial_state,config_text,capacity
from idealAirProperties import IdealAir
from su2IdealAirInterface import MODEL,R,OBSERVED_FLAVOR
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save

KIND='naca-tn3344-spatial-wall-laminar-v1'
FLAVOR='rah-ideal-air-spatial-wall-v1'
SPATIAL=dict(mode='axial-node-temperature-table-v1',schema='RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K',
    environmentVariable='RAH_WALL_TEMPERATURE_FILE',marker='wall',appliedReceipt='rah_wall_temperature_applied.csv',
    interpolation='linear-in-X-no-extrapolation',boundaryFormulaChanged=False,supportsNoTableBaseline=True)
BASE=ROOT/'outputs/v170_validation/mesh02'
PLAN=ROOT/'docs/validation/v171/PLAN.md'
VARIANTS=('observed-scalar','spatial-scalar','spatial-constant','spatial-profile','validation-profile')


def specification(variant):
    assert variant in VARIANTS
    real=variant=='validation-profile'
    return dict(variant=variant,actualValidationMesh=real,updates=2000 if real else 4,every=500 if real else 2,
        maximumSeconds=1800 if real else 30,threads=2,maximumOutputBytes=(1024 if real else 32)*1024**2,
        requiredRuntimeFlavor=OBSERVED_FLAVOR if variant=='observed-scalar' else FLAVOR,
        thermalMode='axialTemperature' if variant in ('spatial-constant','spatial-profile','validation-profile') else 'isothermal',
        nominalScalarTw_K=335.)


def fixture_geometry():
    """模块一：原549点夹具仅合并wall标记，数值组合检查不重造几何。"""
    p,cells,markers,old=small_geometry()
    walls=[face for tag,rows in markers.items() if tag.startswith('wall_') for face in rows]
    markers={tag:rows for tag,rows in markers.items() if not tag.startswith('wall_')};markers['wall']=walls
    labels=np.array(old['nodeIndices'],int)
    ix=np.empty((len(old['offsets']),len(old['s']),len(old['azimuth'])),int)
    for node,(j,i,k) in enumerate(labels):ix[j,i,k]=node
    ix[:,0,:]=ix[:,0,0,None]
    meta={k:old[k] for k in ('s','offsets','azimuth')}
    meta.update(profile='v169-exact-small-geometry-single-wall-v171',readyForCFD=False,newPhysicalQualification=False)
    return p,cells,markers,meta,labels,ix


def spatial_config(variant,c,meta):
    """模块二：数值参数固定，只有有表模式增加官方CustomBoundary标记。"""
    spec=specification(variant)
    text=config_text(c,dict(walls=[dict(marker='wall',givenTw_K=spec['nominalScalarTw_K'])]),spec['updates'],spec['every'],True)
    if spec['thermalMode']=='axialTemperature':text+='MARKER_PYTHON_CUSTOM= (wall)\n'
    return text


def expected_table(variant,meta):
    s=np.array(meta['s'])
    tw=np.full(len(s),335.) if variant=='spatial-constant' else wall_temperature(s)[0]
    x=meridian(s)[:,0]
    return np.column_stack((x[::-1],tw[::-1]))


def source_files():
    """模块三：所有适配/物性/输出/热审核源在启动前固定，任务中不改源。"""
    names=('prepareNACASpatialCase.py','su2NACASpatialInterface.py','prepareSU2NACALaminar.py',
        'naca3344Model.py','su2NACALaminarInterface.py','su2IdealAirInterface.py','su2IdealAirCheckpoints.py',
        'su2ThermalRuntime.py','runEmbeddedSU2.py','nacaSpatialWallFlux.py','analyzeNACASpatialHeat.py')
    return [ROOT/'scripts'/name for name in names]+[PLAN]+[DATA/name for name in
        ('contract.json','source_manifest.json','wall_temperature.csv','heat_parameter.csv')]


def prepare(folder,variant):
    """模块四：独占任务目录、原网格只读硬链接、相容初场及严格具名清单。"""
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs');folder.mkdir(parents=True,exist_ok=False)
    spec=specification(variant);locks={str(p):sha(p) for p in source_files()}
    (folder/'executed_prepare.py').write_bytes(Path(__file__).read_bytes())
    if spec['actualValidationMesh']:
        proof=read(ROOT/'docs/validation/v170/evidence_manifest.json')
        for row in proof['localOnly']:
            p=ROOT/row['path'];assert sha(p)==row['sha256'];locks[str(p)]=row['sha256']
        for name in ('geometry.json','wall_temperature.csv'):
            row=next(r for r in proof['files'] if r['source']=='outputs/v170_validation/mesh02/'+name)
            assert sha(BASE/name)==row['sha256'];locks[str(BASE/name)]=row['sha256']
        from su2NACALaminarInterface import mesh
        p,cells,markers=mesh(BASE/'mesh.su2');meta=read(BASE/'geometry.json')
        for name in ('mesh.su2','mesh_arrays.npz'):os.link(str(BASE/name),str(folder/name))
        with np.load(str(BASE/'mesh_arrays.npz')) as a:labels=a['nodeIndices'];ix=a['columnNodeIDs']
        locks[str(BASE/'geometry.json')]=sha(BASE/'geometry.json')
        assert len(p)==190334 and len(cells)==179568
    else:
        p,cells,markers,meta,labels,ix=fixture_geometry()
        write_mesh(folder/'mesh.su2',p,cells,markers)
        np.savez_compressed(str(folder/'mesh_arrays.npz'),points=p,nodeIndices=labels,columnNodeIDs=ix,
            prisms=np.array([c[1:] for c in cells if c[0]==13]),hexes=np.array([c[1:] for c in cells if c[0]==12]))
    save(folder/'geometry.json',meta)
    c=source_data()[0]['freeStream'];core=IdealAir();locks.update(core.locks)
    initial,thermo=initial_state(p,dict(meta,nodeIndices=labels.tolist()),core,c)
    np.savetxt(str(folder/'initial_flow.csv'),np.column_stack((np.arange(len(p)),p,initial)),delimiter=',',fmt='%.17g',
        comments='',header=','.join('"'+n+'"' for n in ('PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy')))
    (folder/'case.cfg').write_text(spatial_config(variant,c,meta),encoding='ascii')
    profile_sha=None
    if spec['thermalMode']=='axialTemperature':
        if spec['actualValidationMesh']:
            os.link(str(BASE/'wall_temperature.csv'),str(folder/'wall_temperature.csv'))
        else:
            np.savetxt(str(folder/'wall_temperature.csv'),expected_table(variant,meta),delimiter=',',fmt='%.17g',
                comments='',header=SPATIAL['schema'])
        profile_sha=sha(folder/'wall_temperature.csv')
    wallids=np.unique([i for f in markers['wall'] for i in f[1:]])
    contract=dict(schemaVersion=6,initialization=KIND,thermodynamicModel=MODEL,gasConstant=R,
        temperatureRangeK=[60,2000],energyEntropyReference='OTH',turbulent=False,physicsQualified=False,newQualifiedArea=0,
        configurationSHA256=sha(folder/'case.cfg'),meshSHA256=sha(folder/'mesh.su2'),geometrySHA256=sha(folder/'geometry.json'),
        meshArraysSHA256=sha(folder/'mesh_arrays.npz'),checkpointContract=capacity(len(p),len(wallids),len(markers['wall']),
            spec['updates'],spec['every'],spec['maximumOutputBytes']))
    save(folder/'ideal_air_contract.json',contract)
    manifest=dict(caseKind=KIND,specification=spec,requiredRuntimeFlavor=spec['requiredRuntimeFlavor'],
        idealAir=contract,idealAirContractSHA256=sha(folder/'ideal_air_contract.json'),
        thermalBoundary=dict(mode=spec['thermalMode'],givenTw_K=335.,spatialWall=SPATIAL if profile_sha else None),
        thermalProfileSHA256=profile_sha,publicBenchmark=dict(conditions=c,fixture=not spec['actualValidationMesh'],
            experimentalQualification=False,centerAndTransportLimitations='v170 nominal segments; six ambiguous centers; original transport not equivalent'),
        initialization=dict(kind=KIND,priorSuccessfulUpdates=0,initialStateSHA256=sha(folder/'initial_flow.csv'),
            sourceSHA256=locks,initialBlendThickness_m=.001,boundaryLayerSolved=False),
        nodes=len(p),elements=len(cells),expectedVelocity=[-c['speed_m_s'],0,0],
        physicalQualified=False,planSHA256=sha(PLAN))
    save(folder/'case_manifest.json',manifest)
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        maximumSeconds=spec['maximumSeconds'],threads=2,sourceFiles=[p.relative_to(ROOT).as_posix() for p in source_files()])
    if spec['actualValidationMesh']:request['resourceProfile']='v171-naca-spatial-validation'
    save(folder/'request.json',request)
    from su2ThermalRuntime import resolve
    resolve(folder,(folder/'case.cfg').read_text())
    assert all(sha(Path(p))==h for p,h in locks.items())
    print(dict(prepared=True,variant=variant,nodes=len(p),cells=len(cells),outputBudgetBytes=spec['maximumOutputBytes']))
    return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True)
    parser.add_argument('--variant',choices=VARIANTS,required=True);args=parser.parse_args();prepare(args.folder,args.variant)
