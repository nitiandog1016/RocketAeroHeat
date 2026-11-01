"""重放真实隐式增量/松弛，检验只读记录是否改变数值；不把诊断判为物理解。"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from prepareSU2UpdateTraceCase import ROOT,BASE,PLAN,read,sha,config,write
from su2SSTInterface import state

def table(path):
    with path.open() as f:keys=next(csv.reader(f))
    a=np.atleast_2d(np.loadtxt(str(path),delimiter=',',skiprows=1));assert np.isfinite(a).all()
    # RCM本地行序不能作全局身份；必须先验证完整唯一原编号再排序。
    order=np.argsort(a[:,0]);a=a[order]
    assert np.array_equal(a[:,0],np.arange(len(a)))
    return a,keys

def scaled_error(a,b):
    """以各字段最大绝对值归一；幅值小于1时采用绝对误差，避免除零伪放大。"""
    return float(np.max(abs(a-b)/np.maximum(np.max(abs(b),axis=0),1.)))

def replay(folder):
    """按固定CVariable原公式重放；未裁剪量/松弛与最终保存阶段分开核对。"""
    receipt=read(folder/'engine_process.json');assert receipt['status']=='terminal' and receipt['exitCode']==0
    assert all(receipt['inputsStable'].values()) and receipt['physicsModel']['updateTraceStudy']['enabled']
    saved,keys,nd,_,_,_=state(folder/'restart_flow.csv',True,True,True)
    result={};data={};configuration=config(folder/'case.cfg')
    for kind,variables in [('flow',nd+2),('sst',2)]:
        before=folder/(kind+'_update_before.csv');after=folder/(kind+'_update_after.csv')
        for p in [before,after]:assert sha(p)==receipt['outputSHA256'][p.name]
        a,ac=table(before);b,bc=table(after)
        assert len(a)==len(b)==len(saved)
        assert ac==['PointID','DensityUsed','Relaxation']+[prefix+str(v) for prefix in ['Solution_','Old_','Increment_'] for v in range(variables)]
        assert bc==['PointID']+['Solution_'+str(v) for v in range(variables)]
        rho=a[:,1];alpha=a[:,2];current=a[:,3:3+variables];old=a[:,3+variables:3+2*variables];delta=a[:,3+2*variables:]
        assert min(rho)>0 and min(alpha)>=0 and max(alpha)<=1
        if kind=='sst':
            ratio=abs(delta/(rho[:,None]*current+1e-16))
            theoretical=np.minimum(1,float(configuration.get('MAX_UPDATE_SST','1'))/np.maximum(np.max(ratio,axis=1),np.finfo(float).tiny))
            proposed=(rho[:,None]*old+alpha[:,None]*delta)/rho[:,None]
            predicted=np.minimum(np.maximum(proposed,[1e-10,1e-4]),[1e10,1e15])
            output=saved[:,[keys.index('Turb_Kin_Energy'),keys.index('Omega')]]
        else:
            eps=1e-16;rc=np.maximum(current[:,0],eps);rn=np.maximum(rc+delta[:,0],eps)
            ec=current[:,-1]-.5*np.sum(current[:,1:1+nd]**2,axis=1)/rc
            en=current[:,-1]+delta[:,-1]-.5*np.sum((current[:,1:1+nd]+delta[:,1:1+nd])**2,axis=1)/rn
            ratio=np.maximum(abs(delta[:,0])/rc,abs(en-ec)/np.maximum(ec,eps))
            theoretical=np.minimum(1,float(configuration.get('MAX_UPDATE_FLOW','.2'))/np.maximum(ratio,eps))
            proposed=old+alpha[:,None]*delta;predicted=proposed
            output=saved[:,[keys.index('Density')]+[keys.index('Momentum_'+q) for q in 'xyz'[:nd]]+[keys.index('Energy')]]
        theoretical[theoretical<1e-10]=0
        relax_error=float(max(abs(theoretical-alpha)))
        update_error=scaled_error(predicted,b[:,1:]);saved_error=scaled_error(output,b[:,1:])
        assert relax_error<2e-12 and update_error<2e-12
        # 后续流动合法性保护可能回退。明确报告差异，不能强行把两个阶段混为一谈。
        if kind=='sst':assert saved_error<2e-12
        result[kind]=dict(nodes=len(a),relaxationReplayAbsoluteError=relax_error,updateReplayScaledError=update_error,
            savedStateScaledError=saved_error,zeroRelaxationNodes=int(sum(alpha==0)),
            minimumPositiveRelaxation=float(min(alpha[alpha>0])) if any(alpha>0) else None,
            savedStateDifferentNodes=int(sum(np.any(output!=b[:,1:],axis=1))))
        data[kind]=dict(rho=rho,alpha=alpha,current=current,old=old,delta=delta,proposed=proposed,after=b[:,1:],ratios=ratio)
    return result,data,saved,keys

def small():
    target=BASE/'small_comparison_v2.json';assert not target.exists()
    arrays={};locks={};summary={}
    for variant in ['small_official','small_off_v2','small_on_v2']:
        folder=BASE/variant;r=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values()) and audit['interfacePassed']
        a,keys,_,_,_,_=state(folder/'restart_flow.csv',True,True,True);arrays[variant]=a
        assert audit['checkpointLocalIteration']==0
        traces=list(folder.glob('*_update_*.csv'));assert len(traces)==(4 if variant=='small_on_v2' else 0)
        summary[variant]=dict(seconds=r['elapsedSeconds'],exitCode=r['exitCode'],interfacePassed=True)
        for p in [folder/'engine_process.json',folder/'case.cfg',folder/'restart_flow.csv',folder/'sst_state_audit.json',folder/'initial_flow.csv']:
            locks[str(p)]=sha(p)
    assert sha(BASE/'small_official/initial_flow.csv')==sha(BASE/'small_off_v2/initial_flow.csv')==sha(BASE/'small_on_v2/initial_flow.csv')
    assert np.array_equal(arrays['small_off_v2'],arrays['small_on_v2'])
    assert sha(BASE/'small_off_v2/restart_flow.csv')==sha(BASE/'small_on_v2/restart_flow.csv')
    error=scaled_error(arrays['small_off_v2'],arrays['small_official']);assert error<2e-12
    for kind in ['flow','sst']:
        ids=np.loadtxt(str(BASE/'small_on_v2'/(kind+'_update_before.csv')),delimiter=',',skiprows=1,usecols=0)
        assert not np.array_equal(ids,np.arange(len(ids))),'小例必须实际覆盖RCM重编号'
    trace,_,_,_=replay(BASE/'small_on_v2');assert trace['flow']['savedStateScaledError']<2e-12
    result=dict(passed=True,variants=summary,enabledVersusDisabledByteEqual=True,officialVersusCompiledScaledError=error,
        numericalStatementsUnchanged=True,replay=trace,physicsQualified=False,newQualifiedArea=0,
        sourceSHA256=dict(locks,**{str(p):sha(p) for p in [Path(__file__).resolve(),PLAN]}))
    write(target,result);print(json.dumps({k:result[k] for k in ['passed','officialVersusCompiledScaledError','replay']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['small']);p.parse_args();small()
