"""用冻结实际原面数据检查面积身份文案；只更新查看器，保留全部场量和数值字节。"""
from pathlib import Path
import json, shutil, sys
from screenRocketFullFLeastSquares import ROOT, sha


def main(source, target):
    source, target = (Path(p).resolve() for p in (source, target))
    source.relative_to(ROOT.resolve()); target.relative_to(ROOT.resolve())
    assert source.is_dir() and not target.exists()
    shutil.copytree(source, target, ignore=shutil.ignore_patterns(
        'node_wall_velocity_summary.json', 'independent_wall_velocity_view.json'))
    copied = []
    for old in target.rglob('viewer.js'):
        shutil.copy2(ROOT / 'web/viewer.js', old)
        copied.append(str(old.relative_to(target)))
    assert len(copied) == 3
    preserved = []
    for original in (source / 'cases').glob('*/results.json'):
        relative = original.relative_to(source)
        assert original.read_bytes() == (target / relative).read_bytes()
        preserved.append(str(relative))
    receipt = dict(status='built-research-coverage-identity-preview-with-frozen-original-facets',
        sourceFolder=str(source), sourceDataSHA256=sha(source / 'wall-velocity-study-data.json'),
        allOriginalFacetResultsByteIdentical=preserved, updatedViewerAssets=copied,
        viewerSourceSHA256=sha(ROOT / 'web/viewer.js'), newFlowAcceptances=0,
        newActualNSEvaluations=0, newProductionSolvedAreaM2=0, physicalSolution=False)
    (target / 'coverage_identity_source.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf8')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main(*sys.argv[1:3])
