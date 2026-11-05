"""准备 NASA ACC-SBLI 轴对称公开基准，保留原网格和独立模型身份。"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'tests/data/tmr_mach7'
REFERENCE = ROOT / 'tests/data/tmr_mach7_su2'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    Path(path).write_bytes((json.dumps(value, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))


def reference_data():
    """保存官方参考曲线及来源，重用时逐字检查已有清单。"""
    manifest = REFERENCE / 'source_manifest.json'
    if manifest.exists():
        sources = json.loads(manifest.read_text(encoding='utf-8'))
        for item in sources:
            assert sha(REFERENCE / item['path']) == item['sha256']
        return sources
    (REFERENCE / 'raw').mkdir(parents=True, exist_ok=True)
    sources = []
    for part in ('inflow', 'wall'):
        name = 'axiswblim7_wind_' + part + '_sstv.dat'
        url = 'https://tmbwg.github.io/turbmodels/Axiswblim7_validation/' + name
        target = REFERENCE / 'raw' / name
        assert not target.exists(), '已有未登记原始文件，请先核对来源'
        with urllib.request.urlopen(url, timeout=30) as response:
            data = response.read(2 * 1024 * 1024)
        assert data and len(data) < 2 * 1024 * 1024
        assert b'variables' in data.lower() and b'zone' in data.lower()
        target.write_bytes(data)
        sources.append(dict(url=url, path='raw/' + name, bytes=len(data), sha256=sha(target)))
    save_json(manifest, sources)
    return sources


def source_blocks():
    """PLOT3D 的 i 索引最快，单位由厘米转换成米，X 取源 x 的负值。"""
    archive = DATA / 'raw/Axiswblim7_grids.zip'
    with zipfile.ZipFile(str(archive)) as z:
        member = next(n for n in z.namelist() if n.endswith('kussoy_2zn_2levdn.p2dfmt.gz'))
        raw = gzip.decompress(z.read(member))
    numbers = np.fromstring(raw.decode('ascii'), sep=' ')
    assert np.array_equal(numbers[:5], [2, 81, 101, 81, 101])
    n = 81 * 101
    assert len(numbers) == 5 + 4 * n
    blocks = []
    for b in range(2):
        start = 5 + 2 * b * n
        x = numbers[start:start+n].reshape(101, 81)
        r = numbers[start+n:start+2*n].reshape(101, 81)
        blocks.append(np.stack((-x, r), axis=2) * .01)
    return blocks, dict(archive=str(archive.relative_to(ROOT)), archiveSHA256=sha(archive),
                       member=member, decompressedSHA256=hashlib.sha256(raw).hexdigest())


def write_mesh(folder):
    """只合并两块共有接口，并按反射后的正向顺序写出四边形和四类边界。"""
    blocks, source = source_blocks()
    mismatch = float(np.max(np.abs(blocks[0][:, -1] - blocks[1][:, 0])))
    assert mismatch < 2e-14
    points = list(blocks[0].reshape(-1, 2))
    index0 = np.arange(81 * 101).reshape(101, 81)
    index1 = np.empty((101, 81), dtype=int)
    index1[:, 0] = index0[:, -1]
    index1[:, 1:] = np.arange(len(points), len(points) + 80 * 101).reshape(101, 80)
    points.extend(blocks[1][:, 1:].reshape(-1, 2))
    points = np.array(points)
    cells = []
    labels = {}
    for b, ix in enumerate((index0, index1)):
        for j in range(100):
            for i in range(80):
                cells.append([ix[j, i], ix[j+1, i], ix[j+1, i+1], ix[j, i+1]])
        for i in range(80):
            labels[tuple(sorted([ix[0, i], ix[0, i+1]]))] = 'wall'
            labels[tuple(sorted([ix[-1, i], ix[-1, i+1]]))] = 'top'
        for j in range(100):
            col, tag = (0, 'inlet') if b == 0 else (-1, 'outlet')
            labels[tuple(sorted([ix[j, col], ix[j+1, col]]))] = tag
    edges = {}
    counts = Counter()
    for cell in cells:
        for a, b in zip(cell, cell[1:] + cell[:1]):
            key = tuple(sorted((a, b)))
            counts[key] += 1
            edges[key] = (a, b)
    assert set(k for k, n in counts.items() if n == 1) == set(labels)
    boundary = {tag: [] for tag in ('wall', 'inlet', 'top', 'outlet')}
    for key, tag in labels.items():
        boundary[tag].append(edges[key])
    with (folder / 'mesh.su2').open('w', encoding='ascii', newline='\n') as f:
        f.write('NDIME= 2\nNELEM= %d\n' % len(cells))
        for c in cells:
            f.write('9 ' + ' '.join(map(str, c)) + '\n')
        f.write('NPOIN= %d\n' % len(points))
        for i, p in enumerate(points):
            f.write('%.17g %.17g %d\n' % (p[0], p[1], i))
        f.write('NMARK= 4\n')
        for tag, segments in boundary.items():
            f.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n' % (tag, len(segments)))
            for a, b in segments:
                f.write('3 %d %d\n' % (a, b))
    source.update(interfaceMismatchM=mismatch, mergedInterfaceNodes=101,
                  sourceRadiusM=float(blocks[0][0, 0, 1]), nominalRadiusM=.1015)
    return source


def audit_geometry(path):
    """独立读回全部网格，核对局部正性、边配对和旋转体积及壁面面积。"""
    with Path(path).open() as f:
        assert f.readline().strip() == 'NDIME= 2'
        nc = int(f.readline().split('=')[1])
        rows = np.array([list(map(int, f.readline().split())) for _ in range(nc)])
        assert np.all(rows[:, 0] == 9)
        cells = rows[:, 1:]
        npnt = int(f.readline().split('=')[1])
        pts = np.array([list(map(float, f.readline().split())) for _ in range(npnt)])
        assert np.array_equal(pts[:, 2], np.arange(npnt))
        pts = pts[:, :2]
        nm = int(f.readline().split('=')[1])
        boundary = {}
        for _ in range(nm):
            tag = f.readline().split('=')[1].strip()
            n = int(f.readline().split('=')[1])
            e = np.array([list(map(int, f.readline().split())) for _ in range(n)])
            assert np.all(e[:, 0] == 3) and tag not in boundary
            boundary[tag] = e[:, 1:]
        assert not f.read().strip()
    assert nc == 16000 and npnt == 16261 and np.all(pts[:, 1] > 0)
    assert set(boundary) == {'wall', 'inlet', 'top', 'outlet'}
    polygon = pts[cells]
    nextp = np.roll(polygon, -1, axis=1)
    cross = polygon[:, :, 0]*nextp[:, :, 1] - nextp[:, :, 0]*polygon[:, :, 1]
    area = np.sum(cross, axis=1)/2
    turns = np.cross(nextp-polygon, np.roll(polygon, -2, axis=1)-nextp)
    assert np.min(area) > 0 and np.min(turns) > 0
    cell_volume = 2*np.pi*np.sum((polygon[:, :, 1]+nextp[:, :, 1])*cross, axis=1)/6
    assert np.min(cell_volume) > 0
    directed = Counter()
    for c in cells:
        for a, b in zip(c, np.roll(c, -1)):
            directed[(int(a), int(b))] += 1
    assert max(directed.values()) == 1
    exposed = set(e for e in directed if (e[1], e[0]) not in directed)
    boundary_edges = [tuple(e) for seg in boundary.values() for e in seg]
    assert len(boundary_edges) == len(set(boundary_edges))
    assert exposed == set(boundary_edges)
    allseg = np.array(boundary_edges)
    a, b = pts[allseg[:, 0]], pts[allseg[:, 1]]
    boundary_volume = -np.pi*np.sum((b[:, 0]-a[:, 0])*(a[:, 1]**2+a[:, 1]*b[:, 1]+b[:, 1]**2))/3
    volume_error = abs(np.sum(cell_volume)/boundary_volume-1)
    assert volume_error < 2e-10
    wall = pts[boundary['wall']]
    wall_area = float(np.sum(np.pi*np.sum(wall[:, :, 1], axis=1)*np.linalg.norm(wall[:, 1]-wall[:, 0], axis=1)))
    assert np.max(np.ptp(pts[boundary['inlet']][:, :, 0])) < 1e-12
    # 原网格出口带倾斜，完整保留其形状，并检查位于拐角下游。
    assert np.max(pts[boundary['outlet']][:, :, 0]) < 0
    assert np.ptp(pts[boundary['top']][:, :, 1]) < 1e-12
    return dict(passed=True, nodes=npnt, quadrilaterals=nc,
                boundaryCounts={k: len(v) for k, v in boundary.items()},
                minimumPlanarAreaM2=float(min(area)), minimumCornerCrossM2=float(np.min(turns)),
                revolvedVolumeM3=float(boundary_volume), volumeRelativeError=float(volume_error),
                revolvedWallAreaM2=wall_area)


def prepare(folder):
    """固定两更新启动参数，公开基准使用 SST-Vm 独立身份。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT / 'outputs')
    assert not folder.exists(), '保留已有输入和结果，请使用新的计划目录'
    refs = reference_data()
    folder.mkdir(parents=True)
    source = write_mesh(folder)
    geometry = audit_geometry(folder / 'mesh.su2')
    gamma, gas_constant, mach = 1.4, 287.058, 7.11
    total_temperature = 1600*5/9
    total_pressure = 364.8*6894.757293168
    factor = 1 + (gamma-1)*mach**2/2
    temperature = total_temperature/factor
    pressure = total_pressure/factor**(gamma/(gamma-1))
    speed = mach*np.sqrt(gamma*gas_constant*temperature)
    density = pressure/(gas_constant*temperature)
    viscosity = 1.716e-5*(temperature/273.15)**1.5*(273.15+110.4)/(temperature+110.4)
    inlet = '%.17g, %.17g, %.17g, 0, 0' % (temperature, pressure, -speed)
    config = '''SOLVER= RANS
KIND_TURB_MODEL= SST
SST_OPTIONS= (V1994m, VORTICITY)
AXISYMMETRIC= YES
MATH_PROBLEM= DIRECT
RESTART_SOL= NO
MESH_FILENAME= mesh.su2
MESH_FORMAT= SU2
SYSTEM_MEASUREMENTS= SI
REF_DIMENSIONALIZATION= DIMENSIONAL
FLUID_MODEL= IDEAL_GAS
GAMMA_VALUE= 1.4
GAS_CONSTANT= 287.058
INIT_OPTION= TD_CONDITIONS
MACH_NUMBER= 7.11
AOA= 180
SIDESLIP_ANGLE= 0
FREESTREAM_TEMPERATURE= {temperature:.17g}
FREESTREAM_PRESSURE= {pressure:.17g}
FREESTREAM_TURBULENCEINTENSITY= 0.00021
FREESTREAM_TURB2LAMVISCRATIO= 0.001
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716e-5
MU_T_REF= 273.15
SUTHERLAND_CONSTANT= 110.4
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
PRANDTL_LAM= 0.72
TURBULENT_CONDUCTIVITY_MODEL= CONSTANT_PRANDTL_TURB
PRANDTL_TURB= 0.9
MARKER_ISOTHERMAL= (wall, 311)
MARKER_SUPERSONIC_INLET= (inlet, {inlet}, top, {inlet})
MARKER_SUPERSONIC_OUTLET= (outlet)
MARKER_PLOTTING= (wall)
MARKER_MONITORING= (wall)
REF_AREA= {wall_area:.17g}
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CONV_NUM_METHOD_FLOW= AUSM
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG
VENKAT_LIMITER_COEFF= 0.05
TIME_DISCRE_FLOW= EULER_IMPLICIT
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
TIME_DISCRE_TURB= EULER_IMPLICIT
CFL_REDUCTION_TURB= 1.0
LINEAR_SOLVER= BCGSTAB
LINEAR_SOLVER_PREC= LU_SGS
LINEAR_SOLVER_ERROR= 1E-5
LINEAR_SOLVER_ITER= 10
MGLEVEL= 0
CFL_NUMBER= 0.4
CFL_ADAPT= NO
ITER= 2
CONV_FIELD= (RMS_DENSITY, RMS_ENERGY, RMS_TKE, RMS_DISSIPATION)
CONV_RESIDUAL_MINVAL= -8
CONV_STARTITER= 10
TABULAR_FORMAT= CSV
HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF, HEAT, CFL_NUMBER, LINSOL)
CONV_FILENAME= history
SURFACE_FILENAME= surface_flow
RESTART_FILENAME= restart_flow
VOLUME_FILENAME= volume_flow
OUTPUT_FILES= (SURFACE_PARAVIEW_ASCII, RESTART_ASCII)
OUTPUT_WRT_FREQ= 2
WRT_RESTART_OVERWRITE= YES
WRT_SURFACE_OVERWRITE= YES
WRT_RESTART_COMPACT= NO
SCREEN_WRT_FREQ_INNER= 1
VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)
'''.format(temperature=temperature, pressure=pressure, inlet=inlet,
           wall_area=geometry['revolvedWallAreaM2'])
    (folder / 'case.cfg').write_bytes(config.encode('ascii'))
    conditions = dict(mach=mach, gasConstant=gas_constant, gamma=gamma,
                      totalTemperatureK=total_temperature, totalPressurePa=total_pressure,
                      temperatureK=temperature, pressurePa=pressure, densityKgM3=density,
                      speedMS=float(speed), wallTemperatureK=311., turbulenceIntensity=.00021,
                      eddyToLaminarRatio=.001, unitReynoldsPerM=float(density*speed/viscosity),
                      figureUnitReynoldsPerM=5.706e6,
                      unitReynoldsRelativeDifference=float(density*speed/viscosity/5.706e6-1))
    manifest = dict(kind='NASA-ACC-SBLI-M7-axisymmetric-startup-v1',
                    model='SST-Vm (V1994m,VORTICITY)', conditions=conditions,
                    geometry=geometry, source=source, references=refs,
                    boundarySource='NASA NMF and WIND SST-Vm sample',
                    limitations=['two updates only', 'calorically perfect numerical gas',
                                 '20 degree flare at zero incidence', 'no physical qualification'],
                    inputSHA256={n: sha(folder/n) for n in ('case.cfg', 'mesh.su2')},
                    sourceSHA256={str(p.relative_to(ROOT)): sha(p) for p in (
                        Path(__file__).resolve(), DATA/'raw/wind_sstvm.inp', DATA/'raw/sstvm.html',
                        DATA/'raw/kussoy_2zn_4levdn.nmf', REFERENCE/'source_manifest.json')})
    save_json(folder/'benchmark_manifest.json', manifest)
    print(json.dumps(dict(folder=str(folder), geometry=geometry, conditions=conditions)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    prepare(parser.parse_args().folder)
