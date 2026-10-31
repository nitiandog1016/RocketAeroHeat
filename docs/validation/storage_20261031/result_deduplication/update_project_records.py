# -*- coding: utf-8 -*-
"""归档核验完成后更新三份项目记录；本模块不执行数据删除。"""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
FOLDER = Path(__file__).resolve().parent


def main():
    # 模块一：必须已经完成最终核验，并拒绝重复追加执行记录。
    final = json.loads((FOLDER / 'final_postcheck.json').read_text(encoding='utf-8'))
    assert final['status'] == 'passed' and final['originalsRemoved'] == 18
    paths = [ROOT / p for p in ('AGENTS.md', 'docs/OUTPUTS_STORAGE.md', 'docs/PROJECT_GOALS.md')]
    originals = [p.read_bytes() for p in paths]
    marker = '2026-10-31最新存储维护：追加将18份'
    assert marker.encode('utf-8') not in originals[0]
    entry = (
        marker + '已剥离旧矩阵的算子/图结果无损去重至现有mat_history_store_20261031。'
        '原12,132,295,651字节由新增2,481,638,930字节对象和1,163,430字节清单替代，'
        '净减少约9.65GB；执行窗口磁盘空闲实增约9.64GB。18份迁存前MAT的完整SHA均由独立实现'
        '恢复通过后才删除原容器。94项当前依赖（80份MAT）、Q41/Q25、v94—v103及独有结果/'
        '失败/检查点保留；11,371个共享对象通过SHA核对，其他36,673项原文件身份保持，仅用途'
        'README更新且原文另存。全目录最终去重并计入压缩的数据长度约39.92GB，路径累计约'
        '70.02GB不作为实占。38份退役MAT的数值现共用同一对象库，仍只存本地，不能单独删除'
        '共享对象；旧路径旁恢复标记指向清单，按zlib1.2.11恢复到新目录。13项测试及原生删除/'
        '后验通过；见docs/validation/storage_20261031/result_deduplication/。本轮仅存储维护，'
        '不运行MATLAB，不改变物理结果或下一计算起点，不混入OpenCode修改。'
    )
    changes = []
    raw = originals[0]
    nl = b'\r\n' if b'\r\n' in raw else b'\n'
    raw = raw.replace('2026-10-31最新存储维护：20份'.encode('utf-8'),
                      '2026-10-31上一轮存储维护：20份'.encode('utf-8'), 1)
    needle = '## 当前目标与优先级'.encode('utf-8')
    assert raw.count(needle) == 1
    changes.append(raw.replace(needle, entry.encode('utf-8') + nl + nl + needle, 1))

    # 模块二：用途说明注明真实释放量，保留全部历史时点记录。
    raw = originals[1]
    nl = b'\r\n' if b'\r\n' in raw else b'\n'
    heading, rest = raw.split(nl, 1)
    rest = rest.replace('**2026-10-31最新维护已执行并复核：**'.encode('utf-8'),
                        '**2026-10-31上一轮维护已执行并复核：**'.encode('utf-8'), 1)
    lead = (
        '**最新追加维护已执行并复核：** 18份退役算子/图MAT通过独立原SHA恢复校验后，'
        '删除旧容器并共享现有无损归档。原12.13GB仅新增约2.48GB对象和清单，净减少数据长度'
        '约9.65GB；执行窗口磁盘空闲实增约9.64GB。目录去重并计入压缩的数据长度由49.57GB降至'
        '39.92GB。94项当前依赖（80份MAT）、Q41/Q25、v94—v103及全部独有结果/失败/检查点保持；'
        '仅本地用途README另行更新，原文保留。完整说明、各时点实测和恢复命令见'
        '[本轮记录](validation/storage_20261031/result_deduplication/README.md)。38份退役结果'
        '目前共用同一个本地对象库，GitHub没有大型数据备份，不能单独删除共享对象。以下记录'
        '表示各自历史时点。'
    )
    changes.append(heading + nl + nl + lead.encode('utf-8') + nl + rest)

    # 模块三：执行记录明确只做存储维护，科研主线与总体目标继续有效。
    raw = originals[2]
    nl = b'\r\n' if b'\r\n' in raw else b'\n'
    log = (
        '- **2026-10-31用户再次要求清理outputs，追加无损去重已执行并复核：** ' + entry +
        ' 新维护源码、测试、恢复清单和小型回执按约定同步公开GitHub，不上传大型对象库、'
        '不生成工程ZIP或隔离副本。已有区域求解资格保持，下一计算工作仍从v103后的真实20°'
        '横流—收缩—后体、分离/再附及底面工程区域和统一科研UI推进；当前单时刻全面积与UI'
        '最高优先级、完整目标active。'
    )
    changes.append(raw.rstrip(b'\r\n') + nl + nl + log.encode('utf-8') + nl)
    for path, before, after in zip(paths, originals, changes):
        assert path.read_bytes() == before, '文档被并发更新，停止覆盖。'
        path.write_bytes(after)
    print('Updated AGENTS, storage guidance, project execution record')


if __name__ == '__main__':
    main()
