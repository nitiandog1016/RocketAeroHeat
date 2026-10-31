"""独立核对浏览器可见下载链接生成的CSV；不把DOM捕获称为原生下载成功。"""

from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'docs/validation/v104'


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def main():
    # 模块一：读取实际浏览器生成的交线CSV与保存的最终MAT，彼此独立验证。
    csvfile = DEST / 'browser/section_minus1_dom_export.csv'
    rows = list(csv.DictReader(csvfile.open(encoding='utf-8-sig', newline='')))
    matfile = ROOT / 'outputs/v104_validation/actual03/current_case.mat'
    mesh = loadmat(str(matfile), simplify_cells=True)['current']
    vertices = mesh['V']
    triangles = vertices[mesh['F'].astype(int) - 1]
    heat = np.asarray(mesh['heating']['heatFlux']).reshape(-1) / 1000
    selected = np.asarray(mesh['heating']['kinneyAcreage']['acceptedFace']).reshape(-1).astype(bool)
    x = -1.0
    expected = {}
    for face, triangle in enumerate(triangles):
        if not (triangle[:, 0].min() < x < triangle[:, 0].max()):
            continue
        points = []
        for j, k in ((0, 1), (1, 2), (2, 0)):
            a, b = triangle[j], triangle[k]
            if (a[0] - x) * (b[0] - x) < 0:
                points.append(a[1:] + (x - a[0]) / (b[0] - a[0]) * (b[1:] - a[1:]))
        assert len(points) == 2
        expected[face] = points
    assert len(rows) == len(expected) == 92
    assert {int(row['面元编号']) - 1 for row in rows} == set(expected)
    missing = 0
    for row in rows:
        face = int(row['面元编号']) - 1
        assert float(row['截面 X / m']) == x
        a = np.array([float(row['A_Y / m']), float(row['A_Z / m'])])
        b = np.array([float(row['B_Y / m']), float(row['B_Z / m'])])
        points = np.array(expected[face])
        direct = np.max(np.abs(points - np.array([a, b])))
        reverse = np.max(np.abs(points - np.array([b, a])))
        assert min(direct, reverse) < 1e-12
        assert abs(float(row['交线长度 / m']) - np.linalg.norm(a - b)) < 1e-12
        value = row['对流热流密度 / kW/m²']
        if np.isfinite(heat[face]):
            assert abs(float(value) - heat[face]) <= 1e-12 * max(1, abs(heat[face]))
            assert selected[face]
            assert row['状态'] == '工程方法对照'
            assert 'Kinney' in row['算法'] and '入口尚未核验' in row['算法']
            assert '尚未独立验证' in row['适用范围']
        else:
            missing += 1
            assert value == '' and row['状态'] == '未求解'
        assert row['轮廓编号'] == '1' and row['网格环重合时取值侧'] == '+X'
        assert 'M 5' in row['工况'] and 'α 20°' in row['工况']
    assert missing == 34
    report = dict(status='passed-generated-CSV-vs-MAT-and-original-triangle-intersections',
                  captureMethod='read-visible-download-link-data-URL-after-actual-export-button-click',
                  nativeDownloadVerified=False, nativeDownloadAttemptsTimedOut=3,
                  timeoutSecondsPerAttempt=15, reason='IAB-did-not-return-download-event-or-file',
                  rows=len(rows), missingRows=missing, x=x, csvSHA256=sha(csvfile),
                  physicalMAT=matfile.relative_to(ROOT).as_posix(), physicalMATSHA256=sha(matfile),
                  auditorSHA256=sha(Path(__file__)), newPhysicalArea=0)
    target = DEST / 'browser/csv_audit.json'
    raw = (json.dumps(report, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
    if target.exists():
        assert target.read_bytes() == raw
    else:
        target.write_bytes(raw)
    print(json.dumps(dict(status=report['status'], rows=len(rows), missing=missing,
                         nativeDownloadVerified=False)))


if __name__ == '__main__':
    main()
