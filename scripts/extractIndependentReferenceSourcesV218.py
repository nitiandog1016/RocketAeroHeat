"""模块职责：只提取固定构建使用的通量/梯度来源，供完整独立参考的可行性判断。"""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v218_validation'
NAMES=['src/lmr/fluxcalc.d','src/lmr/fvinterface.d','src/lmr/flowstate.d',
       'src/lmr/turbulence.d','src/lmr/flowgradients.d','src/lmr/ufluidblock.d',
       'src/gas/gas_model.d','src/lmr/bc/boundary_interface_effect.d']

def main():
    start=time.monotonic();out=BASE/'source01';out.mkdir(exist_ok=False)
    receipt=ROOT/'outputs/v214_validation/build01/report.json'
    build=json.loads(receipt.read_text());rows=[]
    for rel in NAMES:
        origin='/opt/rocket-aeroheat/eilmer-v5.0.0/source/'+rel
        p=Path(origin);raw=p.read_bytes();h=hashlib.sha256(raw).hexdigest()
        assert h==build['inputSHA256'][origin],origin
        q=out/'source'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(raw)
        rows.append(dict(origin=origin,path=q.relative_to(ROOT).as_posix(),sha256=h,
                         lines=len(raw.splitlines()),bytes=len(raw)))
    report=dict(passed=True,sources=rows,buildReceiptSHA256=hashlib.sha256(receipt.read_bytes()).hexdigest(),
                seconds=time.monotonic()-start,newNativeQueries=0,CFDAdvanceRuns=0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__ == '__main__':
    raise ValueError('v218评估已完成、冻结生产分支关闭；请按NEXT_DECISION登记新任务。')
