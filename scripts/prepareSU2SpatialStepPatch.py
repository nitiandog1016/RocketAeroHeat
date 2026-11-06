"""模块用途：从准确空间Tw源码生成节点限步及查询拒绝上下文，准备过程不覆盖共享树。"""
from pathlib import Path
import hashlib
from buildSU2IdealAirRuntime import ROOT, SOURCE, BUILD, BASE
from prepareSU2ConservativeStepPatch import patch_solver
from su2ThermalRuntime import read, sha

PARENT=ROOT/'outputs/v171_validation/build01'
OLD=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
FLAVOR='rah-ideal-air-spatial-step-v1'
FILE='SU2_CFD/src/variables/CNSVariable.cpp'
SOLVER='SU2_CFD/src/solvers/CNSSolver.cpp'
EULER='SU2_CFD/src/solvers/CEulerSolver.cpp'
ADAPTER='SU2_CFD/include/fluid/CSU2IdealAir.hpp'
TRACE='SU2_CFD/include/variables/ConservativeStepTrace.hpp'
QUERY='SU2_CFD/include/fluid/PropertyQueryContext.hpp'
CHANGED=(FILE,SOLVER,EULER,ADAPTER,TRACE,QUERY)


def replace_once(text,old,new):
    assert text.count(old)==1, '补丁锚点不唯一: '+old[:80]
    return text.replace(old,new,1)


def patch_variable(text):
    """模块一：复用原有限回溯数值规则；只补节点试探/旧态/回溯等拒绝阶段。"""
    text=replace_once(text,'#include "../../include/fluid/CFluidModel.hpp"',
        '#include "../../include/fluid/CFluidModel.hpp"\n#include "../../include/fluid/CSU2IdealAir.hpp"\n'
        '#include "../../include/variables/ConservativeStepTrace.hpp"')
    body=(ROOT/'third_party/ideal_air/ConservativeStep.inc').read_text(encoding='utf8')
    begin=body.index('{')+1
    body=body[:begin]+'''
  // 模块用途：仅定位本次节点查询。局部编号不是全局编号；退出后恢复线程上下文。
  rah_query::Context query_context;
  query_context.phase="node.trial"; query_context.localKnown=true; query_context.local=iPoint;
  rah_query::Guard query_guard(query_context);
'''+body[begin:]
    body=replace_once(body,'    try {\n      query();','    query_guard.phase("node.old");\n    try {\n      query();')
    body=replace_once(body,'    for(unsigned attempt=1;attempt<=32;++attempt) {',
        '    query_guard.phase("node.backtrack");\n    for(unsigned attempt=1;attempt<=32;++attempt) {')
    body=replace_once(body,'    if(!found) {','    if(!found) {\n      query_guard.phase("node.zero-step");')
    body=replace_once(body,'  if(check_dens || check_press || check_sos || check_temp) {',
        '  if(check_dens || check_press || check_sos || check_temp) {\n    query_guard.phase("node.legacy-recovery");')
    start=text.index('bool CNSVariable::SetPrimVar(')
    tail=text.index('  /*--- Set enthalpy ---*/',start)
    return text[:start]+body+text[tail:]


def patch_euler(text):
    """模块二：仅在Riemann物性查询附阶段/原节点身份；不修改特征组合或规定边界压力。"""
    start=text.index('void CEulerSolver::BC_Riemann(')
    end=text.index('\nvoid CEulerSolver::',start+1)
    body=text[start:end]
    marker='    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();'
    body=replace_once(body,marker,marker+'''
    // 模块用途：异常时记录真实边界、局部/全局节点与查询阶段；合法计算不变。
    rah_query::Context query_context;
    query_context.phase="riemann.internal";
    query_context.localKnown=true; query_context.local=iPoint;
    query_context.globalKnown=true; query_context.global=geometry->nodes->GetGlobalIndex(iPoint);
    query_context.markerKnown=true; query_context.marker=val_marker;
    query_context.innerKnown=true; query_context.inner=config->GetInnerIter();
    rah_query::Guard query_guard(query_context);
''')
    marker='      GetFluidModel()->SetTDState_rhoe(Density_i, StaticEnergy_i);'
    body=replace_once(body,marker,marker+'\n      query_guard.phase("riemann.exterior");')
    marker='      GetFluidModel()->SetTDState_rhoe(Density_b, StaticEnergy_b);'
    body=replace_once(body,marker,'      query_guard.phase("riemann.mixed");\n'+marker)
    return text[:start]+body+text[end:]


def prepared_sources():
    """模块三：准确执行源/运行包绑定，未改部分逐字继承，返回覆盖数据而不写共享树。"""
    previous=read(OLD/'runtime.json');tree=read(PARENT/'source_tree.json')
    assert sha(PARENT/'source_tree.json')==previous['sourceTreeSHA256']=='c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8'
    for row in previous['files']:
        assert sha(OLD/row['path'])==row['sha256']
    snapshots={SOURCE/p.relative_to(PARENT/'executed_sources'):p.read_bytes()
        for p in (PARENT/'executed_sources').rglob('*') if p.is_file()}
    for name in (FILE,SOLVER,EULER,ADAPTER):
        p=SOURCE/name
        if p not in snapshots:snapshots[p]=p.read_bytes()
        assert hashlib.sha256(snapshots[p]).hexdigest()==tree[str(p)]
    assert all(hashlib.sha256(value).hexdigest()==tree[str(path)] for path,value in snapshots.items())
    inherited=dict(snapshots)
    snapshots[SOURCE/FILE]=patch_variable(snapshots[SOURCE/FILE].decode('utf8')).encode('utf8')
    snapshots[SOURCE/SOLVER]=patch_solver(snapshots[SOURCE/SOLVER].decode('utf8')).encode('utf8')
    snapshots[SOURCE/EULER]=patch_euler(snapshots[SOURCE/EULER].decode('utf8')).encode('utf8')
    adapter=snapshots[SOURCE/ADAPTER].decode('utf8')
    adapter=replace_once(adapter,'#include "IdealAir.hpp"','#include "IdealAir.hpp"\n#include "PropertyQueryContext.hpp"')
    assert adapter.count('rah_air::state(')==6
    adapter=adapter.replace('rah_air::state(','rah_query::evaluate(')
    snapshots[SOURCE/ADAPTER]=adapter.encode('utf8')
    snapshots[SOURCE/TRACE]=(ROOT/'third_party/ideal_air/ConservativeStepTrace.hpp').read_bytes()
    snapshots[SOURCE/QUERY]=(ROOT/'third_party/ideal_air/PropertyQueryContext.hpp').read_bytes()
    assert {p.relative_to(SOURCE).as_posix() for p in snapshots if snapshots[p]!=inherited.get(p)}==set(CHANGED)
    return previous,tree,inherited,snapshots
