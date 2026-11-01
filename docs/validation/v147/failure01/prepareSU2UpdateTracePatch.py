"""从固定官方归档生成只读更新记录补丁，数值语句逐字保留。"""
from pathlib import Path
import json,hashlib,zipfile,difflib
from prepareSU2CheckpointContinuation import ROOT,sha
OUT=ROOT/'outputs/v147_validation/preflight'
FILES=['SU2_CFD/include/solvers/CScalarSolver.inl','SU2_CFD/include/solvers/CFVMFlowSolverBase.inl']

def block(kind,after):
    flow=kind=='flow'
    condition='R == ENUM_REGIME::COMPRESSIBLE' if flow else 'Conservative'
    opening='if constexpr ('+condition+')' if flow else 'if ('+condition+' && nVar == 2)'
    density='nodes->GetDensity(i)' if flow else 'solver_container[FLOW_SOL]->GetNodes()->GetDensity(i)'
    return '''
  // RAH：全线程同步后只读记录，不写Solution、Jacobian或控制参数。
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    %s {
      RAHWriteUpdateTrace("%s", %s, nPointDomain, nVar, nodes, LinSysSol,
                          [&](unsigned long i) { return %s; });
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
'''%(opening,kind,'true' if after else 'false',density)

def main():
    assert not OUT.exists();OUT.mkdir(parents=True)
    archive=ROOT/'tools/runtime/su2-thermal-build/su2.zip'
    assert sha(archive)=='c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    patch='';records={}
    with zipfile.ZipFile(str(archive)) as z:
        for rel in FILES:
            raw=z.read(next(n for n in z.namelist() if n.endswith('/'+rel)));text=raw.decode('utf-8')
            kind='sst' if 'CScalarSolver' in rel else 'flow'
            start=text.index('::CompleteImplicitIteration(');end=text.index('\n}\n',start)+3
            function=text[start:end]
            anchor='  /*--- Update solution'
            before=function.index(anchor)
            function=function[:before]+block(kind,False)+'\n'+function[before:]
            before=function.index('  for (unsigned short iPeriodic')
            function=function[:before]+block(kind,True)+'\n'+function[before:]
            patched=text[:start]+function+text[end:]
            at=patched.index('\ntemplate <')
            patched=patched[:at]+'#include "RAHUpdateTrace.hpp"\n'+patched[at:]
            # 删除新增记录块与唯一包含后，官方源码必须完整逐字恢复。
            restored=patched.replace('#include "RAHUpdateTrace.hpp"\n','',1)
            for after in [False,True]:restored=restored.replace(block(kind,after)+'\n','',1)
            assert restored==text
            p=OUT/(Path(rel).name+'.patched');p.write_bytes(patched.encode('utf-8'))
            (OUT/Path(rel).name).write_bytes(raw)
            patch+=''.join(difflib.unified_diff(text.splitlines(True),patched.splitlines(True),fromfile='a/'+rel,tofile='b/'+rel))
            records[rel]=dict(originalSHA256=hashlib.sha256(raw).hexdigest(),patchedSHA256=sha(p),numericStatementsUnchanged=True)
    target=ROOT/'third_party/su2_update_trace/update_trace.diff';assert not target.exists();target.write_bytes(patch.encode('utf-8'))
    report=dict(passed=True,files=records,headerSHA256=sha(ROOT/'third_party/su2_update_trace/RAHUpdateTrace.hpp'),patchSHA256=sha(target))
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
