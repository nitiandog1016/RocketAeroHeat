"""模块职责：独立高精度复算已保存局部点云与守恒解码，分离原始量误差；不重跑原生。"""
from pathlib import Path
from decimal import Decimal,localcontext
import csv,json,hashlib,time,bisect
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v217_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_bytes((json.dumps(d,indent=2,allow_nan=False)+'\n').encode('utf8'))
D=lambda v:Decimal.from_float(float(v))

class CaloricReference:
    """用已绑定节点的精确分段Cp积分/二次方程，不调用生产Brent或物性实现。"""
    def __init__(self):
        path=ROOT/'data/thermo/phase_aware_dry_air_v1.json'
        assert sha(path)=='7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'
        v=json.loads(path.read_text());self.t=list(map(D,v['temperatureK']));self.cp=list(map(D,v['cpJkgK']));self.R=D(v['gasConstantJkgK'])
        self.s=[(b-a)/(tb-ta) for a,b,ta,tb in zip(self.cp[:-1],self.cp[1:],self.t[:-1],self.t[1:])]
        h=[Decimal(0)]
        for a,b,ta,tb in zip(self.cp[:-1],self.cp[1:],self.t[:-1],self.t[1:]):h.append(h[-1]+(a+b)*(tb-ta)/2)
        offset=h[self.t.index(D(298.15))];self.e=[v-offset-self.R*t for v,t in zip(h,self.t)]
    def temperature(self,e):
        i=bisect.bisect_right(self.e,e)-1;assert 0<=i<len(self.s)
        cv=self.cp[i]-self.R;delta=e-self.e[i];s=self.s[i]
        dt=2*delta/(cv+(cv*cv+2*s*delta).sqrt())
        assert 0<=dt<=self.t[i+1]-self.t[i]
        return self.t[i]+dt,cv+s*dt
    def primitive(self,U):
        rho=U[0];v=[x/rho for x in U[1:4]];e=U[4]/rho-sum(x*x for x in v)/2
        t,cv=self.temperature(e)
        return [t]+v+[U[5]/rho],cv
    def tangent(self,U,dU):
        q,cv=self.primitive(U);rho=U[0];dr=dU[0]
        dv=[(dU[j+1]-q[j+1]*dr)/rho for j in range(3)]
        de=(dU[4]-U[4]/rho*dr)/rho-sum(q[j+1]*dv[j] for j in range(3))
        return [de/cv]+dv+[(dU[5]-q[4]*dr)/rho]

def main():
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    path=BASE/'control01/stencil.csv';raw=list(csv.DictReader(path.open()))
    rows={}
    for r in raw:
        key=(int(r['regime']),int(r['sample']),int(r['cell']),r['kind'],int(r['index']),int(r['component']))
        assert key not in rows;rows[key]=(float(r['re']),float(r['imScaled']))
    assert all(np.isfinite(v).all() for v in rows.values())
    old=np.genfromtxt(str(ROOT/'outputs/v214_validation/control01/nk_samples.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    baseline={(int(r['regime']),str(r['kind']),int(r['id']),int(r['component'])):r for r in old}
    directions=list(csv.DictReader((ROOT/'outputs/v214_validation/control01/directions.csv').open()))
    dirs={(int(r['regime']),int(r['cell'])):[float(r[k]) for k in ['drho','dmx','dmy','dmz','dE','dsa']] for r in directions}
    discrepancy=0.;nmatch=0
    for key,(value,derivative) in rows.items():
        regime,sample,cell,kind,index,component=key
        if kind=='gradient':
            group='temperatureGradient' if index==0 else 'saGradient' if index==4 else 'velocityGradient'
            cmp=component if index in (0,4) else 3*(index-1)+component;idx=cell
        elif kind=='residual':group='residual';cmp=component;idx=cell
        elif kind=='flux':group='integratedFlux';cmp=component;idx=int(rows[(regime,sample,cell,'faceID',index,0)][0])
        else:continue
        ref=baseline[(regime,group,idx,cmp)]
        a,b=(derivative,float(ref['frozenDerivative'])) if sample==3 else (value,float(ref[['frozen','frozenPlus','frozenMinus'][sample]]))
        discrepancy=max(discrepancy,abs(a-b)/max(1,abs(a),abs(b)));nmatch+=1
    assert discrepancy<=1e-11,discrepancy
    h=D(1e-4);result=[];details=[]
    with localcontext() as ctx:
        ctx.prec=70;gas=CaloricReference()
        for regime in range(3):
            cells=sorted(set(k[2] for k in rows if k[0]==regime));stats=[];decodes=[]
            for cell in cells:
                get=lambda sample,kind,i,j,part=0:D(rows[(regime,sample,cell,kind,i,j)][part])
                weights=[[get(0,'weight',i,j) for i in range(1,6)] for j in range(3)]
                for sample in range(1,4):
                    assert all(get(sample,'weight',i,j)==weights[j][i-1] for i in range(1,6) for j in range(3))
                for field in range(5):
                    for axis in range(3):
                        w=weights[axis]
                        q={s:[get(s,'cloud',i,field) for i in range(6)] for s in range(3)}
                        tangent=[get(3,'cloud',i,field,1) for i in range(6)]
                        apply=lambda a:sum((wi*(a[i+1]-a[0]) for i,wi in enumerate(w)),Decimal(0))
                        gp,gm=apply(q[1]),apply(q[2]);gj=apply(tangent)
                        exactfd=(gp-gm)/(2*h);nativej=get(3,'gradient',field,axis,1)
                        nativefd=(get(1,'gradient',field,axis)-get(2,'gradient',field,axis))/(2*h)
                        scale=max(Decimal(1),abs(nativej),abs(nativefd));gap=abs(nativej-nativefd)/scale
                        # primitive输出量化区间是条件性预算，未包含上游EOS求根/算术误差。
                        radius=lambda q:max(abs(np.nextafter(float(q),np.inf)-float(q)),abs(float(q)-np.nextafter(float(q),-np.inf)))/2
                        quant=Decimal(0)
                        for i,wi in enumerate(w):
                            quant+=abs(wi)*sum((D(radius(q[s][i+1]))+D(radius(q[s][0])) for s in (1,2)),Decimal(0))/(2*h)
                        assembly=max(abs(gp-get(1,'gradient',field,axis)),abs(gm-get(2,'gradient',field,axis)))
                        stat=dict(regime=regime,cell=cell,field=field,axis=axis,originalScaledGap=float(gap),
                                  preciseStencilFDScaledGap=float(abs(nativej-exactfd)/scale),
                                  nativeTangentVsExactWeights=float(abs(nativej-gj)/max(Decimal(1),abs(nativej),abs(gj))),
                                  gradientAssemblyAbsoluteError=float(assembly),
                                  primitiveFDDefectAfterWeights=float(exactfd-gj),
                                  conditionalPrimitiveRoundingRadius=float(quant),
                                  coveredByPrimitiveQuantization=bool(abs(exactfd-gj)<=quant))
                        stats.append(stat)
                U=[get(0,'U0',0,j) for j in range(6)];direction=list(map(D,dirs[(regime,cell)]))
                exactj=gas.tangent(U,direction)
                qp,_=gas.primitive([u+h*d for u,d in zip(U,direction)])
                qm,_=gas.primitive([u-h*d for u,d in zip(U,direction)])
                for field in range(5):
                    nj=get(3,'cloud',0,field,1);nf=(get(1,'cloud',0,field)-get(2,'cloud',0,field))/(2*h)
                    ef=(qp[field]-qm[field])/(2*h)
                    decodes.append(dict(cell=cell,field=field,nativeTangent=float(nj),independentTangent=float(exactj[field]),
                                        tangentScaledDifference=float(abs(nj-exactj[field])/max(Decimal(1),abs(nj),abs(exactj[field]))),
                                        nativeFDMinusTangent=float(nf-nj),preciseFDMinusTangent=float(ef-exactj[field])))
            fail=[s for s in stats if s['originalScaledGap']>1e-5]
            assert max(d['tangentScaledDifference'] for d in decodes)<1e-11
            result.append(dict(regime=regime,cells=cells,gradientComponents=len(stats),originalGradientFailures=len(fail),
                               failuresRemainingAfterExactStencilAssembly=sum(s['preciseStencilFDScaledGap']>1e-5 for s in stats),
                               failuresCoveredByConditionalPrimitiveQuantization=sum(s['coveredByPrimitiveQuantization'] for s in fail),
                               maxTangentVsExactWeights=max(s['nativeTangentVsExactWeights'] for s in stats),
                               maxIndependentPrimitiveTangentDifference=max(d['tangentScaledDifference'] for d in decodes),
                               maxPrimitiveNativeFDDefect=max(abs(d['nativeFDMinusTangent']) for d in decodes),
                               maxPrimitivePreciseFDDefect=max(abs(d['preciseFDMinusTangent']) for d in decodes)))
            details += stats
            save(out/('primitive_regime%d.json'%regime),decodes)
    save(out/'gradients.json',details)
    report=dict(passed=True,observationRows=len(rows),matchedOldComponents=nmatch,selectedOldOutputMaxScaledDifference=discrepancy,
                records=result,independentPrimitiveMapScope='local cell conserved decode with exact binary table Cp integration; not face/boundary/flux or full residual reference',
                primitiveQuantizationBoundScope='conditional on correctly rounded primitive outputs; not a certified upstream evaluation/truncation bound',
                nativePrimitiveTangentProvedLocally=True,fullNumericalAcceptance=False,originalGate=1e-5,perturbationUnchanged=1e-4,
                uniqueWholeResidualCauseProven=False,CFDAdvanceRuns=0,newQualifiedAreaM2=0,
                sourceSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in [path,Path(__file__).resolve(),ROOT/'outputs/v214_validation/control01/nk_samples.csv',ROOT/'data/thermo/phase_aware_dry_air_v1.json']},
                seconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(report))

if __name__ == '__main__':
    raise ValueError('v217观察与审核已结束；当前无构建/查询/网格/CFD准入。')
