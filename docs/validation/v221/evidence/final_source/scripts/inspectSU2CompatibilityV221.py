"""模块职责：定向绑定实际边界与driver接口，避免把共享最新源误作旧对象来源。"""
from pathlib import Path
import hashlib,json,shutil,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v221_validation'
OLD=ROOT/'tools/runtime/su2-thermal-build/source';NEW=ROOT/'outputs/v220_validation/source'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf8')
def main():
    out=BASE/'source01';out.mkdir(parents=True,exist_ok=False);tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json');rows=[]
    names=['SU2_CFD/include/drivers/CDriverBase.hpp','SU2_CFD/include/drivers/CDriver.hpp',
        'SU2_CFD/include/drivers/CSinglezoneDriver.hpp','SU2_CFD/src/drivers/CSinglezoneDriver.cpp',
        'SU2_CFD/src/drivers/CDriver.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp',
        'SU2_CFD/src/solvers/CNSSolver.cpp','SU2_CFD/src/solvers/CTurbSASolver.cpp',
        'SU2_CFD/include/solvers/CScalarSolver.inl','SU2_CFD/include/numerics/turbulent/turb_sources.hpp',
        'SU2_CFD/include/numerics/turbulent/turb_sa_edge_flux.hpp','SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp',
        'SU2_CFD/include/output/CFlowCompOutput.hpp','SU2_CFD/src/output/CFlowCompOutput.cpp']
    for n in names:
        p=NEW/n
        if p.is_file() and str(p.resolve()) in deps:h=deps[str(p.resolve())];assert sha(p)==h
        else:
            h=tree[str(OLD/n)];candidates=[OLD/n]+[ROOT/'outputs'/v/'executed_sources'/n for v in
                ('v167_validation/build01','v164_validation/build01','v163_validation/build03','v165_validation/build01','v154_validation/build01','v158_validation/build01','v162_validation/build01')]
            p=next((q for q in candidates if q.is_file() and sha(q)==h),None);assert p,n
        dest=out/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        rows.append(dict(name=n,path=str(p),sha256=h))
    save(out/'report.json',dict(passed=True,sources=rows,CFDRuns=0))
    print(json.dumps(dict(passed=True,sourceCount=len(rows))))
if __name__=='__main__':main()
