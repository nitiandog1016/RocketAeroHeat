"""显式存储包的准确来源和二进制准入；构建通过不等于数值或物理验证。"""
from pathlib import Path
import json
import re
from runEilmerProcess import COMMIT, sha

PROFILE = 'rah-explicit-dense-v1'
ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'tools/runtime/eilmer-explicit-dense-v1'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def audit_dense_runtime(condition):
    """模块一：禁止使用同名旧安装、修改过的补丁或未正常结束的构建。"""
    d = condition['explicit']
    need(d['profile'] in (PROFILE, 'rah-explicit-dense-development-v1'), '存储包方案错误')
    path = Path(d['runtimePackage']).resolve()
    need(path == PACKAGE / 'runtime.json' and path.is_file(), '存储包位置不符')
    need(sha(path) == d['runtimePackageSHA256'], '存储包清单SHA不符')
    package = json.loads(path.read_text())
    need(package['schema'] == 'RAH.EilmerDenseRuntime.v1' and package['profile'] == PROFILE and
         package['sourceCommit'] == COMMIT and package['boundsCheck'] == 'on', '存储包身份/边界检查不符')
    binary = PACKAGE / 'lmr-run'
    need(Path(package['binary']) == binary, '存储包二进制位置错误')
    expected = {str(PACKAGE / n) for n in ['lmr-run', 'buildinfo.json', 'overlay/src/lmr/fluidblock.d']}
    need(set(package['filesSHA256']) == expected, '存储包文件集合不完整')
    for name, digest in package['filesSHA256'].items():
        need(Path(name).is_file() and sha(name) == digest, '存储包文件改变: ' + name)
    need(sha(package['patch']) == package['patchSHA256'], '存储补丁来源改变')
    receipt = Path(package['buildReceipt'])
    need(sha(receipt) == package['buildReceiptSHA256'], '构建回执改变')
    build = json.loads(receipt.read_text())
    need(build['passed'] and build['nativeExitCode'] == 0 and build['inputsStable'] and
         not build['liveGroupAfterCleanup'], '构建没有正常完成')
    argv = build['argv']
    need('-boundscheck=on' in argv and '-boundscheck=off' not in argv and
         '-of=' + str(binary) in argv and str(PACKAGE / 'overlay/src/lmr/fluidblock.d') in argv,
         '构建命令未采用准确存储实现')
    return package, [path, receipt, Path(package['patch']), *map(Path, package['filesSHA256']), Path(__file__)]


def audit_dense_execution(condition, receipt):
    """模块二：输出不能只凭方案名字认领新实现；核对实际启动路径及锁定SHA。"""
    package, locks = audit_dense_runtime(condition)
    need(len(receipt['commands']) == 1, '存储包运行过程数错误')
    process = receipt['commands'][0]
    need(process['argv'][0] == package['binary'], '实际运行不是登记的存储二进制')
    for path in locks:
        expected = process['inputSHA256'].get(str(path))
        if path == Path(__file__) and expected != sha(path):
            # 模块源码可演化，旧运行必须仍能逐字回到准确执行源；二进制/包不作此例外。
            archived = Path(receipt['case']).parent/'source_at_run'/path.name
            need(archived.is_file() and sha(archived) == expected, '旧存储模块缺少准确执行源')
        else:
            need(expected == sha(path), '存储包运行来源未锁定')
    need(receipt.get('runtimePackageSHA256') == condition['explicit']['runtimePackageSHA256'], '运行缺少存储包身份')
    return True


def audit_dense_allocation(manifest, grid, text):
    """模块三：用原网格独立计数核对实际分配输出，不由补丁文本推定容量已生效。"""
    matches = re.findall(r'^RAH-STORAGE cells=(\d+) ghosts=(\d+) neq=(\d+) nftl=(\d+) scalarBytes=(\d+) stageLength=(\d+) payloadBytes=(\d+)\s*$', text, re.M)
    need(len(matches) == 1, '缺少唯一实际存储容量记录')
    cells, ghosts, neq, levels, width, length, payload = map(int, matches[0])
    need(cells == manifest['mesh']['cells'] and cells == len(grid['cells']), '存储单元计数错误')
    need(ghosts == sum(len(b['ids']) for b in grid['boundaries'].values()), '存储ghost容量不符')
    need((neq, levels, width) == (7, 4, 8), '当前变cp/Wilcox实数四层合同不符')
    need(length == (cells + ghosts)*neq and payload == 2*levels*length*width, '实际数组未按新合同压缩')
    return dict(cells=cells, ghosts=ghosts, equations=neq, timeLevels=levels, scalarBytes=width,
                stageLength=length, arrayPayloadBytes=payload, oldArrayPayloadBytes=payload*levels,
                savedArrayPayloadBytes=payload*(levels-1))
