"""模块职责：只补齐完整组合评审缺少的准确构建源，不运行或修改引擎。"""
from pathlib import Path
import json,time
from inspectSU2RouteV219 import ROOT,BASE,SOURCE,sha,read

NAMES=['SU2_CFD/src/drivers/CDriver.cpp','SU2_CFD/src/integration/CIntegration.cpp',
 'SU2_CFD/include/solvers/CScalarSolver.hpp','SU2_CFD/src/numerics/flow/convection/roe.cpp',
 'SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp','Common/include/CConfig.hpp']

def main():
    start=time.monotonic();out=BASE/'source02';out.mkdir(exist_ok=False)
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');rows=[]
    assert sha(ROOT/'outputs/v167_validation/build01/source_tree.json')=='3b3683a359b1a3f4f5fc6b6f46ac92b5ed3e3c75f5e45c086de6d87e7fd5c540'
    for name in NAMES:
        expected=tree[str(SOURCE/name)]
        candidates=[ROOT/('outputs/'+v+'/executed_sources')/name for v in (
            'v167_validation/build01','v164_validation/build01','v163_validation/build03')]+[SOURCE/name]
        exact=next(p for p in candidates if p.is_file() and sha(p)==expected)
        target=out/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(exact.read_bytes())
        rows.append(dict(relative=name,origin=exact.relative_to(ROOT).as_posix(),sha256=expected,
                         bytes=target.stat().st_size,path=target.relative_to(ROOT).as_posix()))
    report=dict(passed=True,sources=rows,sourceCount=len(rows),seconds=time.monotonic()-start,nativeQueries=0,CFDRuns=0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:report[k] for k in ('passed','sourceCount','seconds')}))

if __name__=='__main__':main()
