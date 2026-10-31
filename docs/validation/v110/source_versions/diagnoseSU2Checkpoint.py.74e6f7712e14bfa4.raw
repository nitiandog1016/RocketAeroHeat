"""核对原生节点残差，并定位层/原壁区/网格质量贡献；不赋予物理资格。"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_volume(path):
    text = Path(path).read_text(encoding='ascii')
    assert '\nASCII\nDATASET UNSTRUCTURED_GRID' in text
    m = re.search(r'POINTS (\d+) \w+\s+(.*?)\s+CELLS', text, re.S)
    n = int(m[1]); xyz = np.fromstring(m[2], sep=' ').reshape(n, 3)
    fields = {}
    pattern = r'(SCALARS|VECTORS) (\S+) \w+(?: (\d+))?\s+(?:LOOKUP_TABLE default\s+)?(.*?)(?=\n(?:SCALARS|VECTORS)|\Z)'
    for m in re.finditer(pattern, text, re.S):
        width = 3 if m[1] == 'VECTORS' else int(m[3])
        value = np.fromstring(m[4], sep=' ')
        assert value.size == n * width, m[2]
        fields[m[2]] = value if width == 1 else value.reshape(n, width)
    assert all(np.all(np.isfinite(v)) for v in fields.values())
    return xyz, fields


def analyze(folder, volume_file):
    receipt = json.loads((folder / 'engine_process.json').read_text())
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0 and all(receipt['inputsStable'].values())
    for name in ('volume_flow.vtk', 'history.csv'):
        assert sha(folder / name) == receipt['outputSHA256'][name]
    v = loadmat(str(volume_file), squeeze_me=True, struct_as_record=False)['volume']
    xyz, fields = read_volume(folder / 'volume_flow.vtk')
    distance, ids = cKDTree(v.vertices).query(xyz)
    assert len(ids) == len(v.vertices) and len(np.unique(ids)) == len(ids) and max(distance) < 1e-6
    inverse = np.argsort(ids); fields = {k: a[inverse] for k, a in fields.items()}
    hist = np.atleast_1d(np.genfromtxt(str(folder / 'history.csv'), delimiter=',', names=True))
    assert len(hist) == 1 and hist['Inner_Iter'][0] == 0
    residual = np.column_stack((fields['Residual_Density'], fields['Residual_Momentum'], fields['Residual_Energy']))
    keys = ['rmsRho', 'rmsRhoU', 'rmsRhoV', 'rmsRhoW', 'rmsRhoE']
    rms = np.log10(np.sqrt(np.mean(residual ** 2, axis=0)))
    error = np.abs(rms - [hist[k][0] for k in keys])
    assert max(error) < 2e-6, '原生空间残差RMS与同次历史不一致'
    wall_ids = np.unique(np.concatenate(v.wall).astype(int)) - 1
    nw = len(wall_ids); assert np.array_equal(wall_ids, np.arange(nw))
    assert len(ids) % nw == 0
    layers = np.arange(len(ids)) // nw
    # 原体网格各层拓扑同构；界面节点按所邻接的不同壁区等权分摊，不任意归入一侧。
    regions = sorted(set(str(x) for x in v.wallRegion)); weights = np.zeros((nw, len(regions)))
    for face, region in zip(v.wall, v.wallRegion):
        weights[np.asarray(face, int) - 1, regions.index(str(region))] = 1
    assert np.all(weights.sum(axis=1) > 0)
    weights /= weights.sum(axis=1)[:, None]
    nodal_weights = np.tile(weights, (len(ids) // nw, 1))
    square = residual ** 2; total = square.sum(axis=0)
    def contributions(mask):
        return dict(zip(keys, (square[mask].sum(axis=0) / total).tolist()))
    region_share = {r: dict(zip(keys, (square * nodal_weights[:, i:i+1]).sum(axis=0) / total)) for i, r in enumerate(regions)}
    layer_share = [dict(layer=int(j), nodes=int(np.sum(layers == j)), residualSquareFraction=contributions(layers == j)) for j in np.unique(layers)]
    near = int(v.meshConstruction.normalLayers)
    top = np.argsort(square[:, -1])[::-1][:30]
    top_nodes = [dict(node=int(i), xyz=v.vertices[i].tolist(), layer=int(layers[i]),
                      wallRegions=[r for j, r in enumerate(regions) if nodal_weights[i, j] > 0],
                      energyResidual=float(residual[i, -1]), energySquareFraction=float(square[i, -1] / total[-1]),
                      orthogonalityDeg=float(fields['Orthogonality'][i]),
                      pressure=float(fields['Pressure'][i]), temperature=float(fields['Temperature'][i]),
                      mach=float(fields['Mach'][i])) for i in top]
    worst = int(np.argmin(fields['Orthogonality']))
    report = dict(status='one-native-sweep-spatial-diagnostic; not-converged-physics',
                  checkpointSourceIteration=3750, diagnosticLocalIteration=0, nodes=len(ids),
                  residualMeaning='SU2 LinSysRes during one new implicit iteration; solution fields are post-update',
                  residualScaleCaveat='integrated control-volume flux imbalance, not volume-normalized error or surface heat-flux error; outer large cells influence this norm',
                  rmsLog10=dict(zip(keys, rms.tolist())), rmsHistoryMaxDifference=float(max(error)),
                  regionAttribution='same meridional topology through all layers; shared region nodes equally partitioned; not surface coverage',
                  regionResidualSquareFraction=region_share, layerResidualSquareFraction=layer_share,
                  wall=contributions(layers == 0), normalShell=contributions((layers > 0) & (layers <= near)),
                  outer=contributions(layers > near), farBoundary=contributions(layers == max(layers)),
                  downstreamOfBody=contributions(v.vertices[:, 0] < np.min(v.sourceVertices[:, 0])),
                  upstreamOfBody=contributions(v.vertices[:, 0] > np.max(v.sourceVertices[:, 0])),
                  lowOrthogonalityBelow15Deg=contributions(fields['Orthogonality'] < 15),
                  minimumOrthogonality=dict(node=worst, valueDeg=float(fields['Orthogonality'][worst]), xyz=v.vertices[worst].tolist(), layer=int(layers[worst])),
                  topEnergyNodes=top_nodes, physicsQualified=False, newQualifiedRocketArea=0,
                  inputSHA256={p.relative_to(ROOT).as_posix(): sha(p) for p in
                               (folder / 'engine_process.json', folder / 'history.csv', folder / 'volume_flow.vtk', volume_file, Path(__file__).resolve())})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True); parser.add_argument('--out', required=True)
    args = parser.parse_args(); dest = ROOT / args.out; assert not dest.exists()
    report = analyze(ROOT / args.folder, ROOT / 'outputs/v108_validation/normal01/normal_volume.mat')
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes((json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf8'))
    print(json.dumps({k: report[k] for k in ['rmsHistoryMaxDifference', 'regionResidualSquareFraction', 'wall', 'normalShell', 'outer', 'farBoundary', 'lowOrthogonalityBelow15Deg', 'minimumOrthogonality']}, indent=2))
    print(json.dumps(report['topEnergyNodes'][:3], indent=2))


if __name__ == '__main__':
    main()
