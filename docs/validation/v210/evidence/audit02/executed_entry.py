"""独立核对原生GasModel的实值、隐函数导数与真实实数差分；不授CFD资格。"""
from pathlib import Path
import json,sys,time,shutil
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir
from phaseAwareDryAir import pressure_ceiling
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v210_validation'
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def inputs():
    """模块一：对既定来流状态增加原生中央差分，不重生成或修改构建输入。"""
    out=BASE/'query_inputs01';out.mkdir(exist_ok=False);prep=BASE/'prepared02'
    config=json.loads((prep/'cases.json').read_text());cases=config['cases']
    a=config['states'][6];p,t=a['p'],a['T'];v=a['values'];rho,u,s,h=v[2],v[3],v[7],v[6]
    pairs={'PT':(p,t),'RU':(rho,u),'RP':(rho,p),'RHOT':(rho,t),'PU':(p,u),'PS':(p,s),'HS':(h,s)}
    scales={'PT':(p,t),'RU':(rho,v[5]*t),'RP':(rho,p),'RHOT':(rho,t),'PU':(p,v[5]*t),'PS':(p,v[4]),'HS':(v[4]*t,v[4])}
    extra=[]
    for mode,(x,y) in pairs.items():
        for axis in (0,1):
            delta=scales[mode][axis]*1e-6
            for sign in (-1,1):
                xx=x+(sign*delta if axis==0 else 0);yy=y+(sign*delta if axis==1 else 0)
                name='fd_%s_%d_%d'%(mode,axis,sign)
                cases.append(dict(id=name,mode=mode,x=xx,y=yy,dx=0.,dy=0.,expected='OK',meta=dict(fd=True,state=6,axis=axis,sign=sign,delta=delta)))
                extra.append('%s %s %.17g %.17g 0 0'%(name,mode,xx,yy))
    for kind in ['real','complex']:(out/(kind+'.tsv')).write_text((prep/(kind+'.tsv')).read_text()+'\n'.join(extra)+'\n')
    save(out/'cases.json',config);shutil.copyfile(__file__,out/'executed_prepare.py')
    save(out/'report.json',dict(passed=True,extraNativeDifferenceRows=28,baseStates=14,CFDRuns=0))
    print('Added 28 native central-difference rows; original preparation unchanged.')

def expected_derivative(model,state,mode,axis):
    """模块二：由独立热力学微分关系建立11量对各入口自变量的解析Jacobian。"""
    p,t=state['p'],state['T'];v=np.array(state['values']);rho,u,cp,cv,h,s,a,mu,k=v[2:]
    i=int(np.searchsorted(model.T,t,side='right')-1);i=min(i,len(model.T)-2)
    dcp=model.slopes[i];dmu=(model.mu[i+1]-model.mu[i])/(model.T[i+1]-model.T[i]);dk=(model.kappa[i+1]-model.kappa[i])/(model.T[i+1]-model.T[i])
    dp=np.array([1,0,rho/p,0,0,0,0,-model.R/p,0,0,0.])
    dt=np.array([0,1,-rho/t,cv,dcp,dcp,cp,cp/t,.5*a*(1/t-model.R*dcp/(cp*cv)),dmu,dk])
    transforms={'PT':[(1,0),(0,1)],'RU':[(model.R*t,0),(rho*model.R/cv,1/cv)],
        'RP':[(0,-t/rho),(1,t/p)],'RHOT':[(model.R*t,0),(rho*model.R,1)],'PU':[(1,0),(0,1/cv)],
        'PS':[(1,model.R*t/(p*cp)),(0,t/cp)],'HS':[(p/(model.R*t),1/cp),(-p/model.R,0)]}
    coefficient=transforms[mode][axis];d=dp*coefficient[0]+dt*coefficient[1]
    characteristic=np.array([p,t,rho,cp*t,cp,cv,cp*t,cp,a,mu,k])
    scales={'PT':(p,t),'RU':(rho,cv*t),'RP':(rho,p),'RHOT':(rho,t),'PU':(p,cv*t),'PS':(p,cp),'HS':(cp*t,cp)}
    scale=np.maximum(np.abs(d),characteristic/scales[mode][axis])
    return d,scale

def audit(branch_recovery=False):
    """模块三：检查完整输出、拒绝类别与误差；任何失败先记录，不覆盖旧回执。"""
    start=time.monotonic();out=BASE/('audit02' if branch_recovery else 'audit01');out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    config=json.loads((BASE/'query_inputs02/cases.json').read_text());model=PhaseAwareViscousAir()
    expected=config['cases'];errors=[];allrows={};maxvalue=0.;maxinverse=0.;maxder=0.;maxfd=0.;derivatives=[]
    counts={};branch_changes=[]
    for variant in ['real','complex']:
        lines=(BASE/'query02'/variant/'stdout.txt').read_text().splitlines()
        if not lines or lines[0]!='MODEL rah-phase-aware-dry-air-eilmer-v1':errors.append([variant,'identity'])
        rows={}
        if variant=='real':
            reused=json.loads((BASE/'query_inputs02/reused_real_rows.json').read_text())
            rows={k:row['line'].split() for k,row in reused.items()}
        for line in lines[1:]:
            fields=line.split()
            if len(fields)<2:errors.append([variant,'malformed',line]);continue
            if fields[0] in rows:errors.append([variant,'duplicate',fields[0]])
            rows[fields[0]]=fields
        subset=[c for c in expected if variant=='complex' or (c['dx']==0 and c['dy']==0)]
        if set(rows)!=set(c['id'] for c in subset):errors.append([variant,'row-set'])
        accepted=0;rejected=0
        for c in subset:
            row=rows.get(c['id'],[])
            if len(row)<2 or row[1]!=c['expected']:
                errors.append([variant,c['id'],'status',row]);continue
            if c['expected']=='REJECT':rejected+=1;continue
            if len(row)!=24:errors.append([variant,c['id'],'column-count',len(row)]);continue
            values=np.array(list(map(float,row[2:]))).reshape(11,2)
            if not np.isfinite(values).all():errors.append([variant,c['id'],'nonfinite']);continue
            rows[c['id']]=values;accepted+=1;meta=c['meta'] or {}
            if meta.get('fd'):continue
            if 'state' in meta:
                state=config['states'][meta['state']];v=np.array(state['values']);t=state['T'];cp=v[4]
                scale=np.abs(v);scale[[3,6]]=cp*t;scale[7]=cp
                value=float(np.max(np.abs(values[:,0]-v)/scale));maxvalue=max(value,maxvalue)
                terr=abs(values[1,0]-t);maxinverse=max(terr,maxinverse)
                if value>2e-10 or terr>2e-7:errors.append([variant,c['id'],'value-or-inverse',value,terr])
                if 'derivativeAxis' in meta:
                    if branch_recovery:
                        # 保留上面的原状态值检查；这里仅选择实际返回状态所在的数学分段。
                        pp,tt=values[0,0],values[1,0]
                        st=model.evaluate(pp,tt);tr=model.interpolated_transport(pp,tt)
                        reference=[pp,tt,float(st['densityKgm3']),float(st['internalEnergyJkg']),float(st['cpJkgK']),float(st['cvJkgK']),
                            float(st['enthalpyJkg']),float(st['entropyJkgK']),float(st['soundSpeedMps']),float(tr['dynamicViscosityPas']),float(tr['conductivityWmK'])]
                        before=int(np.searchsorted(model.T,state['T'],side='right')-1)
                        after=int(np.searchsorted(model.T,tt,side='right')-1)
                        if before!=after:branch_changes.append(dict(id=c['id'],targetTemperatureK=state['T'],returnedTemperatureK=tt,
                            originalInterval=before,actualInterval=after,deltaTemperatureK=tt-state['T']))
                        state=dict(p=pp,T=tt,values=reference)
                    d,scale=expected_derivative(model,state,c['mode'],meta['derivativeAxis'])
                    computed=values[:,1]/1e-30
                    err=float(np.max(abs(computed-d)/scale));maxder=max(maxder,err)
                    derivatives.append(dict(id=c['id'],scaledError=err,expected=d.tolist(),observed=computed.tolist()))
                    if err>2e-8:errors.append([variant,c['id'],'derivative',err])
            else:
                if c['id']=='minimum':
                    if values[1,0]!=60. or values[3,0]>=0 or values[0,0]>pressure_ceiling(60.):errors.append([variant,'minimum'])
                elif c['id'].startswith('limit'):
                    if values[0,0]!=float(c['x']) or values[1,0]!=float(c['y']):errors.append([variant,c['id'],'changed-state'])
        counts[variant]=dict(expected=len(subset),accepted=accepted,rejected=rejected)
        allrows[variant]=rows
    # 用原生实数helper的正负扰动差分，独立检查原生复数结果。
    for c in expected:
        meta=c['meta'] or {}
        if not meta.get('fd') or meta['sign']!=1:continue
        mode=c['mode'];axis=meta['axis'];plus=allrows['real'].get(c['id']);minus=allrows['real'].get('fd_%s_%d_-1'%(mode,axis))
        if not isinstance(plus,np.ndarray) or not isinstance(minus,np.ndarray):continue
        fd=(plus[:,0]-minus[:,0])/(2*meta['delta']);d,scale=expected_derivative(model,config['states'][6],mode,axis)
        err=float(np.max(abs(fd-d)/scale));maxfd=max(maxfd,err)
        if err>2e-5:errors.append(['finite-difference',mode,axis,err])
    report=dict(passed=not errors,errors=errors,counts=counts,maximumScaledValueError=maxvalue,
        maximumTemperatureInverseErrorK=maxinverse,maximumScaledComplexDerivativeError=maxder,
        maximumScaledNativeDifferenceError=maxfd,derivativeCases=len(derivatives),
        elapsedSeconds=time.monotonic()-start,nativeFactoryDispatched=True,originalQueryPlanPassed=False,
        composedFromOriginalAndExplicitRecovery=True,engineBuilds=0,CFDRuns=0,
        physicalQualification=False,derivatives=derivatives,branchChanges=branch_changes,
        derivativeReference='actual verified returned state with right interval at exact knots' if branch_recovery else 'original intended state')
    save(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['derivatives','errors','branchChanges']}))
    print(json.dumps(errors[:24]));return not errors

if __name__=='__main__':
    if sys.argv[1:]==['inputs']:inputs()
    elif sys.argv[1:]==['audit']:raise SystemExit(0 if audit() else 1)
    elif sys.argv[1:]==['audit-branches']:raise SystemExit(0 if audit(True) else 1)
    else:raise ValueError('inputs or audit')
