"""显式工作区包准入：绑定实际构建、真实长度和单元对象，禁止旧包冒认。"""
from pathlib import Path
import json
import re
from runEilmerProcess import COMMIT, sha
from eilmerDenseRuntime import audit_dense_allocation

PROFILE = 'rah-explicit-workspace-v1'
ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT/'tools/runtime/eilmer-explicit-workspace-v1'


def need(ok, message):
    if not ok: raise ValueError(message)


def audit_workspace_runtime(condition):
    """模块一：完整独立包及准确构建来源；不关闭改变湍流处理的NK编译版本。"""
    d=condition['explicit'];need(d['profile']==PROFILE,'工作区方案不符')
    path=Path(d['runtimePackage']).resolve()
    need(path==PACKAGE/'runtime.json' and path.is_file(),'工作区包位置错误')
    need(sha(path)==d['runtimePackageSHA256'],'工作区包SHA错误')
    package=json.loads(path.read_text())
    need(package['schema']=='RAH.EilmerWorkspaceRuntime.v1' and package['profile']==PROFILE and
        package['sourceCommit']==COMMIT and package['boundsCheck']=='on','工作区包身份错误')
    binary=PACKAGE/'lmr-run';need(Path(package['binary'])==binary,'工作区二进制位置不符')
    expected={str(PACKAGE/n) for n in ['lmr-run','buildinfo.json','overlay/src/lmr/fluidblock.d','overlay/src/lmr/ufluidblock.d']}
    need(set(package['filesSHA256'])==expected,'工作区文件集合不符')
    for name,digest in package['filesSHA256'].items():
        need(Path(name).is_file() and sha(name)==digest,'工作区文件改变: '+name)
    need(sha(package['patch'])==package['patchSHA256'],'工作区补丁改变')
    receipt=Path(package['buildReceipt']);need(sha(receipt)==package['buildReceiptSHA256'],'工作区构建回执改变')
    build=json.loads(receipt.read_text())
    need(build['passed'] and build['nativeExitCode']==0 and build['inputsStable'] and not build['liveGroupAfterCleanup'],'工作区构建失败')
    argv=build['argv']
    need('-d-version=newton_krylov' in argv and '-boundscheck=on' in argv and '-boundscheck=off' not in argv and
        '-of='+str(binary) in argv and all(str(PACKAGE/'overlay/src/lmr'/n) in argv for n in ['fluidblock.d','ufluidblock.d']),
        '工作区实际编译参数不符')
    return package,[path,receipt,Path(package['patch']),*map(Path,package['filesSHA256']),Path(__file__),Path(__file__).with_name('eilmerDenseRuntime.py')]


def audit_workspace_execution(condition, receipt):
    """模块二：来源由实际命令和锁定SHA证明，不由显示标签认领。"""
    package,locks=audit_workspace_runtime(condition)
    need(len(receipt['commands'])==1,'工作区运行过程数错误')
    process=receipt['commands'][0]
    need(process['argv'][0]==package['binary'],'未使用工作区二进制')
    for path in locks:
        need(process['inputSHA256'].get(str(path))==sha(path),'工作区运行来源未锁定')
    need(receipt.get('runtimePackageSHA256')==condition['explicit']['runtimePackageSHA256'],'运行包身份缺失')
    return True


def audit_workspace_allocation(manifest,grid,text):
    """模块三：所有备份长度须为0，正常梯度/时间层不删，已构造单元无隐式对象。"""
    storage=audit_dense_allocation(manifest,grid,text)
    matches=re.findall(r'^RAH-WORKSPACE jacobian=(\d+) sourceBackup=(\d+) viscousBackup=(\d+) interpBackup=(\d+) perturbBackup=(\d+) viscous=(\d+) interp=(\d+) interpWorkspace=(\d+)$',text,re.M)
    need(len(matches)==1,'缺少唯一工作区实际长度记录')
    values=list(map(int,matches[0]));n=storage['cells']+storage['ghosts']
    need(values[:5]==[0]*5 and values[5:]==[n]*3,'工作区未省去备份或误删正常梯度')
    cells=re.findall(r'^RAH-WORKSPACE-CELLS internal=(\d+) implicitObjects=(\d+)$',text,re.M)
    need(len(cells)==1 and tuple(map(int,cells[0]))==(storage['cells'],0),'实际单元仍含隐式对象')
    return dict(**storage,workspaceLengths=values,internalCellsChecked=storage['cells'],implicitObjects=0,
        dynamicScope='one registered explicit small case; not all configurations or whole-domain memory fit')
