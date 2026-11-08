"""模块职责：用既存结果分离输出量化、散度组装与点云条件数；不重跑原生算子。"""
from pathlib import Path
from decimal import Decimal,localcontext
import json,time,csv
import numpy as np
from readEilmerNative import read_grid,sha
from auditFrozenMappedGeometryV213 import FACES
from eilmerNKSAResourceBudget import ROOT,save
BASE=ROOT/'outputs/v216_validation';OLD=ROOT/'outputs/v214_validation'

def main():
    start=time.monotonic();out=BASE/'numerical01';out.mkdir(exist_ok=False)
    data=OLD/'control01/nk_samples.csv';raw=np.genfromtxt(str(data),delimiter=',',names=True,dtype=None,encoding='utf8')
    grid=read_grid(OLD/'control01/lmrsim/grid/grid-0000')
    geometry=np.genfromtxt(str(ROOT/'outputs/v213_validation/control03/geometry.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    cells=geometry[geometry['kind']=='cell'];faces=geometry[geometry['kind']=='face'];vol=cells['volumeOrArea']
    xyz=lambda x:np.column_stack([x[k] for k in ['x','y','z']])
    lookup={tuple(sorted(row[row>=0])):i for i,row in enumerate(grid['faces'])}
    ids=np.array([[lookup[tuple(sorted(row[list(local)]))] for local in FACES] for row in grid['cells'][:,:6]])
    dx=xyz(faces)[ids]-xyz(cells)[:,None,:]
    signs=np.sign(np.einsum('ijk,ijk->ij',dx,np.column_stack([faces[k] for k in ['nx','ny','nz']])[ids]))
    assert np.all(abs(signs)==1)
    # 模块一：准确p=.75加权矩阵。SVD只分析几何敏感性，不替换生产梯度。
    theta=np.linalg.norm(dx,axis=2)**(-.75);A=dx*theta[:,:,None]
    sv=np.linalg.svd(A,compute_uv=False);condition=sv[:,0]/sv[:,-1]
    weights=np.linalg.pinv(A)*theta[:,None,:]
    full_weights=np.concatenate([-weights.sum(axis=2)[:,:,None],weights],axis=2)
    gain=np.sum(abs(full_weights),axis=2)
    # gamma_7约束五项已保存积分通量的累加、倒体积和乘法；不含上游通量/物性评估误差。
    u=2.0**-53;gamma=7*u/(1-7*u);h=1e-4;records=[]
    with (out/'components.csv').open('w',newline='') as fp:
        writer=csv.writer(fp);writer.writerow(['regime','kind','cell','component','J','FD','scaledGap','outputQuantizationRadius','assemblyAllowance'])
        for regime in range(3):
            r=raw[raw['regime']==regime]
            result=dict(regime=regime,groups={})
            for kind in ['residual','temperatureGradient','velocityGradient','saGradient','integratedFlux']:
                a=r[r['kind']==kind];plus=a['frozenPlus'];minus=a['frozenMinus'];j=a['frozenDerivative']
                fd=(plus-minus)/(2*h);scale=np.maximum(1,np.maximum(abs(fd),abs(j)))
                # 下一/上一浮点间隔取较大者，覆盖2的幂处的不对称间隔。
                radius=lambda q:np.maximum(np.nextafter(q,np.inf)-q,q-np.nextafter(q,-np.inf))/2
                quant=(radius(plus)+radius(minus))/(2*h)+u*abs(fd)*2
                err=abs(j-fd)/scale;failed=err>1e-5
                result['groups'][kind]=dict(count=len(a),failed=int(failed.sum()),maxScaledGap=float(err.max()),
                    finalOutputQuantizationRadiusMax=float((quant/scale).max()),
                    failuresCoveredByFinalOutputRounding=int(np.sum(failed&(abs(j-fd)<=quant))),
                    residualFailureNotRelabeled=True)
                allowance=np.zeros(len(a))
                if kind=='residual':
                    flux=r[r['kind']=='integratedFlux'];pp=flux['frozenPlus'].reshape(-1,6);mm=flux['frozenMinus'].reshape(-1,6)
                    # 前五式源=0；SA源单列，不把线性散度的界套入非线性源项。
                    boundp=gamma*np.sum(abs(pp[ids]),axis=1)/vol[:,None]
                    boundm=gamma*np.sum(abs(mm[ids]),axis=1)/vol[:,None]
                    allowance=((boundp+boundm)/(2*h)).reshape(-1)
                    exactp=np.zeros((len(cells),5));exactm=exactp.copy()
                    with localcontext() as ctx:
                        ctx.prec=70
                        for ci in range(len(cells)):
                            for eq in range(5):
                                for values,dest in [(pp,exactp),(mm,exactm)]:
                                    s=sum((Decimal.from_float(float(values[fi,eq]))*Decimal(int(-sgn)) for fi,sgn in zip(ids[ci],signs[ci])),Decimal(0))
                                    dest[ci,eq]=float(s/Decimal.from_float(float(vol[ci])))
                    savedp=plus.reshape(-1,6)[:,:5];savedm=minus.reshape(-1,6)[:,:5]
                    bp=boundp[:,:5];bm=boundm[:,:5]
                    # 高精度结果转回double也有最后一次量化，单列加入条件性组装界。
                    delta=np.maximum(abs(savedp-exactp)/(bp+radius(exactp)),abs(savedm-exactm)/(bm+radius(exactm)))
                    fd_from_flux=(exactp-exactm)/(2*h);jj=j.reshape(-1,6)[:,:5]
                    resdiff=abs(jj-fd_from_flux)/np.maximum(1,np.maximum(abs(jj),abs(fd_from_flux)))
                    result['assembly']=dict(firstFiveSavedStateCoveredByGamma7=bool(np.all(delta<=1)),
                        maxFractionOfConditionalBound=float(delta.max()),
                        finiteDifferenceGapAfter70DigitAssembly=float(resdiff.max()),
                        firstFiveFailedAfter70DigitAssembly=int(np.sum(resdiff>1e-5)),
                        boundScope='already-rounded integrated fluxes; no flux/gradient/EOS error or FD truncation bound',
                        SAIncluded=False)
                    # 最坏误差指定到实际原cell/母cell，供下一次最小可区分观察直接使用。
                    worst=int(np.argmax(err));result['worstResidual']=dict(cell=int(a['id'][worst]),equation=int(a['component'][worst]),
                        scaledGap=float(err[worst]),J=float(j[worst]),FD=float(fd[worst]),
                        volumeM3=float(vol[int(a['id'][worst])]),faceIDs=ids[int(a['id'][worst])].tolist())
                for k in np.flatnonzero(failed):writer.writerow([regime,kind,int(a['id'][k]),int(a['component'][k]),j[k],fd[k],err[k],quant[k],allowance[k]])
            records.append(result)
    report=dict(passed=True,schema='RAH.FrozenDifferenceResolution.v1',records=records,
        weightedQRPointCloudCondition=dict(min=float(condition.min()),max=float(condition.max()),median=float(np.median(condition)),
            largestComponentInputGainPerM=float(gain.max()),smallestSingularValue=float(sv[:,-1].min()),fullColumnRank=True),
        originalFiniteDifferenceAcceptance=False,fullNumericalAcceptance=False,uniqueCauseProven=False,
        missingEvidence=['certified error of pre-assembly real flux/gradient/EOS evaluations',
            'complete FD truncation bound or independent full directional reference'],
        originalGate=1e-5,perturbationUnchanged=h,nativeRuns=0,newQualifiedAreaM2=0,
        sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__).resolve(),data,
            ROOT/'outputs/v213_validation/control03/geometry.csv',
            ROOT/'outputs/v211_validation/source01/source/lmr/flowgradients.d',
            ROOT/'tools/runtime/eilmer-nk-geometric-v1/overlay/src/lmr/fluidfvcell.d']},seconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
