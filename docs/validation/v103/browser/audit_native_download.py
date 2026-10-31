"""核对真正下载到磁盘的36条周向记录及可追溯数值，不依赖DOM预览换行。"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent


def main():
    data_file = ROOT / 'outputs/v103_validation/view01/workbench/ale-data.json'
    data = json.loads(data_file.read_text(encoding='utf-8'))
    case, frame = data['cases'][2], data['cases'][2]['frames'][1]
    chosen = [i for i, kind in enumerate(case['boundaryKind']) if kind == 2]
    vertices = np.asarray(frame['position'])
    faces = np.asarray(case['faces'], int)[chosen]
    file = HERE / 'native_section.csv'
    rows = list(csv.reader(file.read_text(encoding='utf-8-sig').splitlines()))
    assert len(rows) == 37 and len(rows[0]) == 18
    observed = set()
    for row in rows[1:]:
        assert len(row) == 18 and float(row[0]) == .04
        local = int(row[2]) - 1
        assert 0 <= local < len(chosen) and local not in observed
        observed.add(local)
        triangle = vertices[faces[local]]
        for point in ([.04, float(row[5]), float(row[6])], [.04, float(row[7]), float(row[8])]):
            coefficient = np.linalg.lstsq((triangle[1:] - triangle[0]).T, np.asarray(point) - triangle[0], rcond=None)[0]
            assert coefficient.min() > -1e-9 and coefficient.sum() < 1 + 1e-9
            assert np.linalg.norm(triangle[0] + coefficient @ (triangle[1:] - triangle[0]) - point) < 2e-12
        assert float(row[10]) == frame['wallPressurePa'][chosen[local]]
        assert row[12] == case['methods']['wallPressurePa']
        assert case['source']['evidenceSHA256'] in row[14]
        assert float(row[14].split('t=')[1].split(' s')[0]) == frame['timeS']
    report = dict(status='passed-actual-browser-download-barycentric-fields-and-current-time', rows=36,
                  xM=.04, timeS=frame['timeS'], fileSHA256=hashlib.sha256(file.read_bytes()).hexdigest(),
                  dataSHA256=hashlib.sha256(data_file.read_bytes()).hexdigest(),
                  productionQualified=False)
    (HERE / 'native_download_audit.json').write_bytes((json.dumps(report, indent=2) + '\n').encode('utf-8'))
    print(report)


if __name__ == '__main__':
    main()
