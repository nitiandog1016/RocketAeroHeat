"""真实入口启动前准入拒绝检查；只读已有输入，不运行CFD。"""
from pathlib import Path
import copy,json,sys,time,shutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runOpenFoamRocketDevelopment import BASE,validate_preparation,read_json,control_text
from readOpenFoamRocketState import mapped_geometry


def main():
    start=time.monotonic();out=BASE/'admission_tests01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_test.py')
    budget=read_json(BASE/'budget01/report.json');review=read_json(ROOT/'outputs/v200_validation/rocket_relocation_audit01/report.json')
    assert validate_preparation(budget,review)
    mapped_geometry(ROOT/'outputs/v200_validation/rocket03')
    rows=[]
    edits=[('failed_budget',lambda b,r:b.update(passed=False)),('wrong_cells',lambda b,r:b.update(cellCount=14400)),
        ('missing_original_wall',lambda b,r:b.update(originalWallFaces=240)),('changed_area',lambda b,r:b.update(originalAreaM2=2)),
        ('changed_Co',lambda b,r:b.update(targetMaxCo=.3)),('wrong_dt',lambda b,r:b.update(initialTargetDeltaTSeconds=1e-7)),
        ('failed_geometry',lambda b,r:r.update(passed=False)),('boundary_changed',lambda b,r:r.update(originalWallExact=False)),
        ('connections_changed',lambda b,r:r.update(connectivityExact=False)),('missing_sha',lambda b,r:b.update(sourceSHA256={})),
        ('wrong_source_sha',lambda b,r:b['sourceSHA256'].update({next(iter(b['sourceSHA256'])):'0'*64}))]
    for name,edit in edits:
        b=copy.deepcopy(budget);r=copy.deepcopy(review);edit(b,r)
        try:validate_preparation(b,r)
        except ValueError as e:rows.append(dict(name=name,rejected=True,reason=str(e)));continue
        raise AssertionError('未拒绝 '+name)
    c=control_text()
    for token in ['endTime 0.00045;','maxCo 0.2;','writeInterval 0.000225;','startTime 1e-8; endTime 2e-8;']:
        assert token in c
    result=dict(passed=True,positiveChecks=3,rejectionCount=len(rows),rejections=rows,CFDRuns=0,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['passed','positiveChecks','rejectionCount','elapsedSeconds']}))


if __name__=='__main__':main()
