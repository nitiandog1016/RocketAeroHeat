"""复用原3750检查点，准备一次有预算的一层粗网格对照，不覆写历史输入。"""
from pathlib import Path
import hashlib
import json
import os

ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    source=ROOT/'outputs/v109_validation/rocket_seed01'
    control=ROOT/'outputs/v111_validation/control01'
    out=ROOT/'outputs/v112_validation/mg01';assert not out.exists()
    old=json.loads((source/'engine_process.json').read_text(encoding='utf-8'))
    assert old['status']=='terminal' and all(old['inputsStable'].values())
    assert sha(source/'restart_flow.csv')==old['outputSHA256']['restart_flow.csv']
    cfg=(control/'case.cfg').read_text(encoding='ascii');assert cfg.count('MGLEVEL= 0')==1
    cfg=cfg.replace('MGLEVEL= 0','MGLEVEL= 1')
    cfg+='MGCYCLE= V_CYCLE\nMG_PRE_SMOOTH= (1, 2)\nMG_POST_SMOOTH= (0, 0)\nMG_CORRECTION_SMOOTH= (0, 0)\nMG_DAMP_RESTRICTION= 0.5\nMG_DAMP_PROLONGATION= 0.5\n'
    out.mkdir(parents=True)
    for oldname,newname in [('mesh.su2','mesh.su2'),('restart_flow.csv','initial_flow.csv')]:
        os.link(str(source/oldname),str(out/newname))
    (out/'case.cfg').write_bytes(cfg.encode('ascii'))
    inputs=[source/'restart_flow.csv',source/'mesh.su2',source/'engine_process.json',control/'case.cfg',
            ROOT/'docs/validation/v112/PLAN.md',Path(__file__).resolve()]
    manifest=dict(sourceCheckpointIteration=3750,control='outputs/v111_validation/control01',
                  method='one-coarse-level FAS V cycle; unchanged finest-grid physical model and original surface',
                  qualificationThresholdUnchanged=-8,requiredLogImprovement=.3,physicsQualified=False,
                  inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs})
    (out/'study_manifest.json').write_bytes((json.dumps(manifest,indent=2)+'\n').encode('utf-8'))
    req=dict(workDirectory=str(out),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=360,
             sourceFiles=['scripts/prepareSU2MultigridStudy.py','docs/validation/v112/PLAN.md',
                          'outputs/v112_validation/mg01/study_manifest.json'])
    (out/'engine_request.json').write_bytes((json.dumps(req,indent=2)+'\n').encode('utf-8'))
    print(str(out/'engine_request.json'))


if __name__=='__main__':main()
