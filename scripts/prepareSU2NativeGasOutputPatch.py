"""模块用途：准确v171父源的纯输出补丁；准备只返回字节，不覆盖共享源码。"""
import hashlib
from pathlib import Path

from buildSU2IdealAirRuntime import ROOT, SOURCE, BUILD, BASE
from nacaConstantCpModel import MODEL, GAMMA, R
from su2ThermalRuntime import read, sha

PARENT = ROOT/'outputs/v171_validation/build01'
OLD = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
FLAVOR = 'rah-native-ideal-gas-spatial-wall-v1'
FILE = 'SU2_CFD/src/output/CFlowOutput.cpp'
OBJECT = 'SU2_CFD/src/libSU2core.a.p/output_CFlowOutput.cpp.obj'
TREE_SHA = 'c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8'
RUNTIME_SHA = 'e794985c282f73913c679eb72ad1818c9d56e97371813270ec7364865bfa3b3d'
PARENT_OBJECTS_SHA = '56904a4c83263d5eaed66219b9d425c52b5425f3211c010ba29a8b0851fd0ecb'
PARENT_COMMANDS_SHA = 'a896ed2ecd5283d5d7ee29147c912f6270d12cec567960b63c87908351465e04'
CACHE_CANDIDATES = (
    'Common/src/libSU2Common.a.p/geometry_CGeometry.cpp.obj',
    'SU2_CFD/src/libSU2core.a.p/python_wrapper_structure.cpp.obj',
    'SU2_CFD/src/libSU2core.a.p/solvers_CTurbSSTSolver.cpp.obj',
    'SU2_CFD/src/libSU2core.a.p/variables_CNSVariable.cpp.obj',
)
OUTPUT_CONTRACT = dict(
    schema='RAH.NativeIdealGas.FirstCompletedOutput.v1',
    supportedModel='IDEAL_GAS', thermodynamicModel=MODEL, gamma=GAMMA, gasConstant=R,
    supportedCaseKind='naca-tn3344-constant-cp-spatial-wall-v1',
    steadyFirstCompletedUpdateLocalIteration=0,
    iterationMeaning='local-zero-is-first-completed-update-not-initial-state',
    changedSource=FILE, changedObject=OBJECT, baseCOutputChanged=False,
    firstFrameOutputModelExtension='IDEAL_GAS',
    existingFrequencyAndForcedFinalRulesUnchanged=True,
    numericalUpdateSourceUnchanged=True, sameNumericalUpdateSources=True,
    outputScheduleChanged=True, nativeSharedFrameEqualityVerified=False,
    selectedCaseNumericalEqualityRequiresSeparateRunEvidence=True,
    asciiSignificantDigits=17, physicalQualified=False,
)


def patch_output(text):
    """模块一：仅扩大首帧输出条件，不改基类、输出字段及任何求解更新。"""
    old = 'return ((Iter > 0 || config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP) &&'
    new = ('return ((Iter > 0 || config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP ||\n'
           '             config->GetKind_FluidModel() == IDEAL_GAS) &&')
    assert text.count(old) == 1, '准确首帧补丁锚点不唯一'
    assert text.index(old) > text.index('bool CFlowOutput::WriteVolumeOutput(')
    answer = text.replace(old, new, 1)
    assert answer.replace(new, old, 1) == text
    return answer


def prepared_sources():
    """模块二：绑定准确父包与19份执行源，返回唯一改变的覆盖数据。"""
    assert sha(OLD/'runtime.json') == RUNTIME_SHA
    previous = read(OLD/'runtime.json')
    assert previous['flavor'] == 'rah-ideal-air-spatial-wall-v1'
    assert previous['asciiSignificantDigits'] == 17
    assert sha(PARENT/'source_tree.json') == previous['sourceTreeSHA256'] == TREE_SHA
    assert sha(PARENT/'objects_after.json') == PARENT_OBJECTS_SHA
    assert sha(PARENT/'commands.json') == PARENT_COMMANDS_SHA
    for item in previous['files']:
        path = (OLD/item['path']).resolve(); path.relative_to(OLD.resolve())
        assert sha(path) == item['sha256'], '准确父运行包发生变化'
    tree = read(PARENT/'source_tree.json')
    inherited = {SOURCE/path.relative_to(PARENT/'executed_sources'): path.read_bytes()
        for path in (PARENT/'executed_sources').rglob('*') if path.is_file()}
    assert len(inherited) == 19 and SOURCE/FILE in inherited
    assert all(hashlib.sha256(raw).hexdigest() == tree[str(path)] for path, raw in inherited.items())
    snapshots = dict(inherited)
    snapshots[SOURCE/FILE] = patch_output(inherited[SOURCE/FILE].decode('utf8')).encode('utf8')
    assert {str(p.relative_to(SOURCE)).replace('\\','/') for p in snapshots
            if snapshots[p] != inherited[p]} == {FILE}
    return previous, tree, inherited, snapshots


def object_plan():
    """模块三：逐项核对298对象；非父对象只能从已确认的父SHA缓存恢复。"""
    assert sha(PARENT/'objects_after.json') == PARENT_OBJECTS_SHA
    parent = read(PARENT/'objects_after.json')
    before = {str(path.relative_to(BUILD)): sha(path) for path in BUILD.rglob('*.obj')}
    assert len(parent) == len(before) == 298 and set(parent) == set(before)
    available = {str(Path(name)): PARENT/'produced_objects'/name for name in CACHE_CANDIDATES}
    assert all(sha(path) == parent[name] for name, path in available.items()), '父对象缓存身份不符'
    different = {name for name in parent if before[name] != parent[name]}
    assert different.issubset(set(available) | {str(Path(OBJECT))}), '未登记的非父对象；拒绝混包'
    cached = {name: available[name] for name in different if name != str(Path(OBJECT))}
    return parent, before, cached
