"""模块用途：独立准备NACA半球柱三维层流入口；不修改原箭体或旧轴对称合同。"""
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np

from analyzeSU2DualVolumes import FACES
from idealAirProperties import IdealAir
from meshSU2HybridVolume import wedge_jacobian_min
from naca3344Model import DATA, D, R, SHOULDER, END, meridian, source_data, wall_temperature, strip_mean
from su2IdealAirCheckpoints import output_plan
from su2IdealAirInterface import AXISYMMETRIC_FLAVOR, OBSERVED_FLAVOR, OBSERVATION, MODEL, R as GAS_R
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

KIND = 'naca-tn3344-quarter-3d-laminar-v1'
STAGE = 'same-flow-update-laminar-conservative-and-primitive-v1'


def small_geometry():
    """模块一：四分之一真实鼻尖，用棱柱接六面体；不含1/r轴对称源。"""
    s = np.r_[np.linspace(0, SHOULDER, 9), np.linspace(SHOULDER, END, 5)[1:]]
    y = np.r_[0., 2e-6*5.**np.arange(8)]
    phi = np.linspace(0, np.pi/2, 5)
    points, labels = [], []
    ix = np.empty((len(y), len(s), len(phi)), int)
    for j, offset in enumerate(y):
        for i, length in enumerate(s):
            xr = meridian(length, offset)
            for k, angle in enumerate(phi):
                if i == 0 and k:
                    ix[j,i,k] = ix[j,i,0]
                    continue
                ix[j,i,k] = len(points)
                points.append([xr[0], xr[1]*np.cos(angle), xr[1]*np.sin(angle)])
                labels.append([j,i,k])
    p = np.array(points)
    p[np.abs(p)<1e-17] = 0
    cells = []
    for j in range(len(y)-1):
        for i in range(len(s)-1):
            for k in range(len(phi)-1):
                if i == 0:
                    ids = [ix[j,0,0],ix[j,1,k],ix[j,1,k+1],
                           ix[j+1,0,0],ix[j+1,1,k],ix[j+1,1,k+1]]
                    cells.append([13]+ids)
                else:
                    ids = [ix[j,i,k],ix[j,i+1,k],ix[j,i+1,k+1],ix[j,i,k+1],
                           ix[j+1,i,k],ix[j+1,i+1,k],ix[j+1,i+1,k+1],ix[j+1,i,k+1]]
                    cells.append([12]+ids)
    lookup = defaultdict(list)
    for cell in cells:
        for local in FACES[cell[0]]:
            face = [cell[1+i] for i in local]
            lookup[tuple(sorted(face))].append(face)
    markers = defaultdict(list)
    labels = np.array(labels)
    for variants in lookup.values():
        assert len(variants) <= 2
        if len(variants) == 2:
            a,b = variants
            ar = a[a.index(min(a)):]+a[:a.index(min(a))]
            br = b[b.index(min(b)):]+b[:b.index(min(b))]
            assert ar == [br[0]]+br[:0:-1]
            continue
        face = variants[0]
        jl,il,kl = labels[face].T
        if np.all(jl == 0):
            tag = 'wall_%03d'%min(il)
        elif np.all(jl == len(y)-1):
            tag = 'outer'
        elif np.all(il == len(s)-1):
            tag = 'outlet'
        elif np.all(np.abs(p[face,2]) < 1e-15):
            tag = 'sym_z'
        elif np.all(np.abs(p[face,1]) < 1e-15):
            tag = 'sym_y'
        else:
            raise AssertionError('Unclassified complete boundary')
        markers[tag].append([5 if len(face)==3 else 9]+face)
    walls = []
    for i in range(len(s)-1):
        Tw,area = strip_mean(lambda a: wall_temperature(a)[0], s[i], s[i+1])
        walls.append(dict(marker='wall_%03d'%i, lowerS_m=float(s[i]), upperS_m=float(s[i+1]),
            givenTw_K=Tw, exactFullRevolvedArea_m2=area,
            extension=bool(s[i+1] > source_data()[1]['s_over_D'][-1]*D)))
    return p, cells, dict(sorted(markers.items())), dict(s=s.tolist(), offsets=y.tolist(),
        azimuth=phi.tolist(), nodeIndices=labels.tolist(), walls=walls)


def geometry_audit(p, cells, markers, meta):
    """模块二：完整边界/子体积、棱柱全参考正性和解析回转域对照。"""
    volumes = []
    area = 0.
    minimum_sub = float('inf')
    for cell in cells:
        q = p[cell[1:]]
        center = q.mean(axis=0)
        vol = 0.
        for local in FACES[cell[0]]:
            face = q[local]
            fc = face.mean(axis=0)
            for a,b in zip(face,np.roll(face,-1,axis=0)):
                v = np.dot(a-center,np.cross(b-center,fc-center))/6
                minimum_sub = min(minimum_sub,float(v))
                vol += v
        assert vol > 0 and minimum_sub > 0
        volumes.append(vol)
    w = np.array([c[1:] for c in cells if c[0] == 13])
    prism_jac = float(min(wedge_jacobian_min(p[w])))
    assert prism_jac > 0
    bvol = 0.
    for name, rows in markers.items():
        for row in rows:
            face = p[row[1:]]
            fc = face.mean(axis=0)
            for a,b in zip(face,np.roll(face,-1,axis=0)):
                normal = np.cross(b-a,fc-a)/2
                bvol += np.dot(a,normal)/3
                if name.startswith('wall_'):
                    area += np.linalg.norm(normal)
    closure = abs(bvol/sum(volumes)-1)
    assert closure < 2e-10
    height = meta['offsets'][-1]
    smooth = (2*np.pi/3*((R+height)**3-R**3)+np.pi*((R+height)**2-R**2)*(.127-R))/4
    ideal_area = (2*np.pi*R**2+2*np.pi*R*(.127-R))/4
    assert np.array_equal(p[0], [0,0,0])
    return dict(passed=True,nodes=len(p),cells=len(cells),prisms=len(w),hexes=len(cells)-len(w),
        wallFaces=sum(len(v) for k,v in markers.items() if k.startswith('wall_')),
        minimumMedianSubVolume_m3=minimum_sub, minimumPrismJacobian=prism_jac,
        fluidVolume_m3=float(sum(volumes)), boundaryVolumeRelativeError=float(closure),
        smoothGeometryVolumeRelativeError=float(abs(sum(volumes)/smooth-1)),
        wallArea_m2=area, smoothWallArea_m2=ideal_area, quarterSectorRadians=np.pi/2,
        containsActualNose=True, axisymmetricSourceEnabled=False, physicalQualified=False)


def write_mesh(path,p,cells,markers):
    """模块三：保持所有物理marker和原节点编号，独占新目录。"""
    with path.open('w',encoding='ascii',newline='\n') as f:
        f.write('NDIME= 3\nNELEM= %d\n'%len(cells))
        for c in cells:f.write(' '.join(map(str,c))+'\n')
        f.write('NPOIN= %d\n'%len(p))
        for i,pt in enumerate(p):f.write('%.17g %.17g %.17g %d\n'%(*pt,i))
        f.write('NMARK= %d\n'%len(markers))
        for name,faces in markers.items():
            f.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n'%(name,len(faces)))
            for face in faces:f.write(' '.join(map(str,face))+'\n')


def initial_state(p,meta,core,c):
    """模块四：无SST的相容正状态；固定C1壁距混合只是初猜，不是已发展边界层。"""
    labels=np.array(meta['nodeIndices'])
    s=np.array(meta['s'])[labels[:,1]]
    distance=np.array(meta['offsets'])[labels[:,0]]
    tw,_=wall_temperature(s)
    eta=np.minimum(distance/.001,1.)
    blend=eta**2*(3-2*eta)
    wall=core.state('PT',c['pressure_Pa'],tw)
    target=(1-blend)*wall['h']+blend*c['totalEnthalpy_OTH_Jkg']
    u=np.column_stack((-c['speed_m_s']*blend,np.zeros((len(p),2))))
    thermo=core.state('Ph',c['pressure_Pa'],target-.5*np.sum(u*u,axis=1))
    rho=thermo['rho']
    conservative=np.column_stack((rho,rho[:,None]*u,rho*(thermo['e']+.5*np.sum(u*u,axis=1))))
    return conservative,thermo


def config_text(c,meta,iterations,every,high_precision=False):
    walls=meta['walls']
    tags=', '.join(w['marker'] for w in walls)
    temperatures=', '.join('%s, %.17g'%(w['marker'],w['givenTw_K']) for w in walls)
    # 此入口只开放已登记的小检查规模；完整候选必须另行计划和网格误差评估。
    text=f'''SOLVER= NAVIER_STOKES
KIND_TURB_MODEL= NONE
AXISYMMETRIC= NO
MATH_PROBLEM= DIRECT
RESTART_SOL= YES
MESH_FILENAME= mesh.su2
MESH_FORMAT= SU2
SYSTEM_MEASUREMENTS= SI
REF_DIMENSIONALIZATION= DIMENSIONAL
FLUID_MODEL= IDEAL_AIR_VARIABLE_CP
GAS_CONSTANT= {GAS_R:.17g}
INIT_OPTION= TD_CONDITIONS
MACH_NUMBER= {c['Mach']:.17g}
AOA= 180
SIDESLIP_ANGLE= 0
FREESTREAM_TEMPERATURE= {c['staticTemperature_K']:.17g}
FREESTREAM_PRESSURE= {c['pressure_Pa']:.17g}
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716e-5
MU_T_REF= 273.15
SUTHERLAND_CONSTANT= 110.4
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
PRANDTL_LAM= 0.72
MARKER_ISOTHERMAL= ({temperatures})
MARKER_SUPERSONIC_INLET= (outer, {c['staticTemperature_K']:.17g}, {c['pressure_Pa']:.17g}, {-c['speed_m_s']:.17g}, 0, 0)
MARKER_SUPERSONIC_OUTLET= (outlet)
MARKER_SYM= (sym_y, sym_z)
MARKER_PLOTTING= ({tags})
MARKER_MONITORING= ({tags})
REF_AREA= {np.pi*R*R/4:.17g}
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG
VENKAT_LIMITER_COEFF= 0.05
TIME_DISCRE_FLOW= EULER_IMPLICIT
LINEAR_SOLVER= BCGSTAB
LINEAR_SOLVER_PREC= LU_SGS
LINEAR_SOLVER_ERROR= 1E-5
LINEAR_SOLVER_ITER= 10
MGLEVEL= 0
CFL_NUMBER= 0.4
CFL_ADAPT= NO
ITER= {iterations}
CONV_FIELD= (RMS_DENSITY, RMS_ENERGY)
CONV_RESIDUAL_MINVAL= -8
CONV_STARTITER= 10
TABULAR_FORMAT= CSV
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF, HEAT, CFL_NUMBER, LINSOL)
CONV_FILENAME= history
SURFACE_FILENAME= surface_flow
RESTART_FILENAME= restart_flow
VOLUME_FILENAME= volume_flow
OUTPUT_FILES= (SURFACE_PARAVIEW_ASCII, RESTART_ASCII)
OUTPUT_WRT_FREQ= {every}
WRT_RESTART_OVERWRITE= NO
WRT_SURFACE_OVERWRITE= NO
WRT_RESTART_COMPACT= NO
SCREEN_WRT_FREQ_INNER= 100
VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL)
SOLUTION_FILENAME= initial_flow
READ_BINARY_RESTART= NO
'''
    if high_precision:
        # 已有17位包的成熟一般EOS特征边界，单独声明，不冒称与旧外推算子相同。
        lines=[line for line in text.splitlines() if not line.startswith(
            ('MARKER_SUPERSONIC_INLET=', 'MARKER_SUPERSONIC_OUTLET='))]
        lines.append('MARKER_RIEMANN= (outer, STATIC_SUPERSONIC_INFLOW_PT, %.17g, %.17g, -1, 0, 0, outlet, STATIC_PRESSURE, %.17g, 0, 0, 0, 0)'%
            (c['pressure_Pa'],c['staticTemperature_K'],c['pressure_Pa']))
        text='\n'.join(lines)+'\n'
    return text


def capacity(nodes,wall_nodes,wall_faces,iterations,every,maximum_bytes):
    value=output_plan(nodes,wall_nodes,wall_faces,iterations,every,maximum_bytes)
    value['stage']=STAGE
    return value


def prepare(folder):
    """模块五：写具名小例合同，不改默认，不启动计算，也不绕过公共监督器。"""
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    folder.mkdir(parents=True,exist_ok=False)
    p,cells,markers,meta=small_geometry()
    audit=geometry_audit(p,cells,markers,meta)
    dataset,_,_=source_data();c=dataset['freeStream']
    core=IdealAir()
    initial,thermo=initial_state(p,meta,core,c)
    write_mesh(folder/'mesh.su2',p,cells,markers)
    headers=['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy']
    np.savetxt(str(folder/'initial_flow.csv'),np.column_stack((np.arange(len(p)),p,initial)),
        delimiter=',',fmt='%.17g',comments='',header=','.join('"'+n+'"' for n in headers))
    (folder/'case.cfg').write_text(config_text(c,meta,3,1,True),encoding='ascii')
    save(folder/'geometry.json',meta)
    save(folder/'geometry_audit.json',audit)
    wallids=np.unique([i for k,v in markers.items() if k.startswith('wall_') for row in v for i in row[1:]])
    locks=dict(core.locks)
    for file in (Path(__file__).resolve(),ROOT/'scripts/naca3344Model.py',DATA/'contract.json',
                 DATA/'source_manifest.json',DATA/'wall_temperature.csv',DATA/'heat_parameter.csv'):
        locks[str(file)]=sha(file)
    init=dict(kind=KIND,priorSuccessfulUpdates=0,sourceSHA256=locks,
        initialStateSHA256=sha(folder/'initial_flow.csv'),temperatureRangeK=[float(min(thermo['T'])),float(max(thermo['T']))],
        boundaryLayerSolved=False,initialBlendThickness_m=.001)
    contract=dict(schemaVersion=5,initialization=KIND,thermodynamicModel=MODEL,
        gasConstant=GAS_R,temperatureRangeK=[60,2000],energyEntropyReference='OTH',
        turbulent=False,physicsQualified=False,newQualifiedArea=0,
        configurationSHA256=sha(folder/'case.cfg'),meshSHA256=sha(folder/'mesh.su2'),
        geometrySHA256=sha(folder/'geometry.json'),
        checkpointContract=capacity(len(p),len(wallids),audit['wallFaces'],3,1,32*1024**2))
    save(folder/'ideal_air_contract.json',contract)
    manifest=dict(caseKind=KIND,requiredRuntimeFlavor=OBSERVED_FLAVOR,runtimeObservation=OBSERVATION,
        runtimeUse='existing-general-3d-NS-with-axisymmetric-source-disabled',
        idealAir=contract,idealAirContractSHA256=sha(folder/'ideal_air_contract.json'),initialization=init,
        nodes=len(p),elements=len(cells),publicBenchmark=dict(fixture=True,conditions=c,
            newPhysicalQualification=False,comparisonReady=False,normalization='local-Nu-over-sqrt-Re-not-dimensional-experiment'),
        expectedVelocity=[-c['speed_m_s'],0,0],physicalQualified=False)
    save(folder/'case_manifest.json',manifest)
    save(folder/'request.json',dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        maximumSeconds=30,threads=2,sourceFiles=['scripts/prepareSU2NACALaminar.py','scripts/naca3344Model.py',
            'scripts/su2NACALaminarInterface.py','docs/validation/v169/PLAN.md']))
    from su2ThermalRuntime import resolve
    resolve(folder,(folder/'case.cfg').read_text())
    return audit


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',required=True)
    args=parser.parse_args()
    print(prepare(args.folder))
