"""模块职责：只用已保存primitive重建基准守恒量，先检验完整参考的输入是否可精确恢复。"""
from pathlib import Path
import csv,json,hashlib,struct,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v218_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_bytes((json.dumps(d,indent=2,allow_nan=False)+'\n').encode('utf8'))

def main():
    start=time.monotonic();out=BASE/'reconstruct01';out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'outputs/v214_validation/control01/directions.csv'
    data=list(csv.DictReader(source.open()));observed=ROOT/'outputs/v217_validation/control01/stencil.csv'
    actual={}
    for row in csv.DictReader(observed.open()):
        if row['sample']=='0' and row['kind']=='U0':actual[(int(row['regime']),int(row['cell']),int(row['component']))]=float(row['re'])
    restored=[];mismatch=[];checked=0
    for row in data:
        regime,cell=int(row['regime']),int(row['cell']);r,u,x,y,z,nu=[float(row[k]) for k in ('rho','u','vx','vy','vz','nuhat')]
        # 模块一：严格沿固定encode的标量顺序，不引入FMA/重排或改变能量零点。SA不含k。
        ke=0.5*((x*x+y*y)+z*z)
        U=[r,r*x,r*y,r*z,r*(u+ke),r*nu]
        restored.append([regime,cell]+U)
        for j,v in enumerate(U):
            key=(regime,cell,j)
            if key in actual:
                checked+=1;bits=lambda x:struct.pack('>d',x).hex()
                if bits(v)!=bits(actual[key]):mismatch.append(dict(regime=regime,cell=cell,component=j,reconstructed=v,observed=actual[key],reconstructedBits=bits(v),observedBits=bits(actual[key])))
    assert len(data)==1248 and checked==96
    with (out/'conserved.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['regime','cell','rho','mx','my','mz','E','rhoSA']);w.writerows(restored)
    report=dict(passed=not mismatch,checkedComponents=checked,checkedStates=16,fullStatesReconstructed=len(data),
                bitwiseMismatches=mismatch,allOtherStatesIndividuallyObserved=False,
                noNewNativeQueries=True,fullNumericalAcceptance=False,
                sourceSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in [source,observed,Path(__file__).resolve(),ROOT/'tools/runtime/eilmer-nk-geometric-v1/overlay/src/lmr/fluidfvcell.d']},
                seconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(report));return not mismatch

if __name__=='__main__':raise SystemExit(0 if main() else 1)
