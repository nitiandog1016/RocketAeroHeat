"""实数算子尺度的独立解析真值、旧配置兼容与拒绝；不是原生流动/实验验证。"""
from pathlib import Path
import argparse,copy,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareEilmerRocketCase import validate_condition,make_job
from readEilmerNative import audit_numerics


def test(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    old=ROOT/'outputs/v180_validation/small01'
    m=json.loads((old/'case_manifest.json').read_text());c=copy.deepcopy(m['condition'])
    assert make_job(m['mesh']['nativeSourcePath'],c)==(old/'job.lua').read_text()
    c['development']['profile']='rah-nk-real-v1';validate_condition(c)
    job=make_job(m['mesh']['nativeSourcePath'],c)
    assert job.split('NewtonKrylovGlobalConfig')[0]==(old/'job.lua').read_text().split('NewtonKrylovGlobalConfig')[0]
    nk=json.loads((old/'lmrsim/nk-config').read_text());a=2.0**-26
    nk.update(use_scaling=True,use_real_valued_frechet_derivative=True,frechet_derivative_perturbation=a,preconditioner_perturbation=a)
    audit_numerics(nk,c);rejected=[]
    for key,value in [('use_scaling',False),('use_real_valued_frechet_derivative',False),
                      ('frechet_derivative_perturbation',1e-30),('preconditioner_perturbation',1e-30),
                      ('use_preconditioner',False),('preconditioner','jacobi')]:
        bad=copy.deepcopy(nk);bad[key]=value
        try:audit_numerics(bad,c)
        except ValueError:rejected.append(key)
        else:raise AssertionError(key)
    # 模块：量级跨越的守恒量风格向量，解析二次残差给出独立Jacobian。
    # 使用固定方向，不扫描参数；只校验差分误差与数值可表示性，不声称真实NS算子。
    scale=np.array([.02,30,10,1,20000,.04,200.]);u=scale*np.array([1.,1.,1.,0.,1.,1.,1.])
    z=scale*np.array([.4,-.3,.2,.1,-.5,.6,-.4])
    def residual(v):
        w=v/scale
        return w*w+.3*np.sum(w)
    exact=(np.diag(2*u/scale)+.3*np.ones((7,7)))/scale[np.newaxis,:]
    sigma=np.mean(a*(np.abs(u)+1))/max(np.linalg.norm(z),np.sqrt(1e-7))
    jv=(residual(u+sigma*z)-residual(u))/sigma
    error=float(np.linalg.norm(jv-exact@z)/np.linalg.norm(exact@z))
    columns=[]
    for i in range(7):
        delta=a*(abs(u[i])+1);v=u.copy();v[i]+=delta
        assert v[i]!=u[i]
        columns.append((residual(v)-residual(u))/delta)
    matrix=np.stack(columns,axis=1)
    matrix_error=float(np.linalg.norm((matrix-exact)*scale)/np.linalg.norm(exact*scale))
    assert error<1e-5 and matrix_error<1e-5
    assert np.array_equal(u+1e-30*z,u)==False  # 零分量可能改变，不能说所有分量都不变。
    erased=int(np.count_nonzero((u+1e-30*z)==u));assert erased>=6
    (out/'condition.json').write_text(json.dumps(c,indent=2)+'\n')
    (out/'job.lua').write_text(job)
    r=dict(passed=True,legacyDevelopmentTextIdentical=True,physicalLuaPrefixIdentical=True,
           independentAnalyticJvRelativeError=error,independentScaledJacobianRelativeError=matrix_error,
           realPerturbation=a,complexSizeErasedComponents=erased,rejected=rejected,nativeRuns=0,
           scope='Independent analytic scale check, not native Jacobian or physical qualification',
           sourceSHA256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
                         ['tests/testEilmerRealProfile.py','scripts/prepareEilmerRocketCase.py','scripts/readEilmerNative.py']},
           elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');return r


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);print(json.dumps(test(p.parse_args().out)))
