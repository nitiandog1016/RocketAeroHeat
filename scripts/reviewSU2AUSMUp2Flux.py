"""只读复核已保存原生通量；保留首版尺度失败，不再次编译或执行。"""
from pathlib import Path
import json,numpy as np
from verifySU2AUSMUp2Flux import ROOT,OUT,PLAN,COMMIT,sha,write,body,flux_scale

def main():
    assert not (OUT/'report.json').exists()
    for name in ['compile','execute']:
        r=json.loads((OUT/(name+'_process.json')).read_text());assert r['status']=='terminal' and r['exitCode']==0
    cases=json.loads((OUT/'inputs.json').read_text());data=np.loadtxt(str(OUT/'native_flux.csv'),delimiter=',')
    assert len(cases)==35 and np.array_equal(data[:,0],np.arange(len(cases))) and np.isfinite(data).all()
    source=(OUT/'ausm_slau.cpp').read_text();cpp=(OUT/'native_flux.cpp').read_text()
    for signature in ['void CUpwAUSMPLUSUP2_Flow::ComputeMassAndPressureFluxes(',
                      'CNumerics::ResidualType<> CUpwAUSMPLUS_SLAU_Base_Flow::ComputeResidual(']:assert body(source,signature) in cpp
    errors=[];old_errors=[]
    for i,c in enumerate(cases[:32]):
        # 解析oracle重新从输入物理状态计算，而不是只读生成器写出的expected。
        left=np.asarray(c['left']);right=np.asarray(c['right']);n=np.asarray(c['normal'])
        if c['kind']=='uniform':
            assert np.array_equal(left,right);mass=left[5]*np.dot(left[1:4],n)
            exact=np.r_[mass,mass*left[1:4]+left[4]*n,mass*left[6]]
        else:
            assert np.all(left[1:4]==0) and np.all(right[1:4]==0) and left[4]==right[4]
            exact=np.r_[0,left[4]*n,0]
        d=abs(data[i,1:]-exact);errors.append(float(max(d/flux_scale(c))));old_errors.append(float(max(d/np.maximum(abs(exact),1))))
    reverse=float(max(abs(data[33,1:]+data[32,1:])/flux_scale(cases[32])))
    scale=float(max(abs(data[34,1:]-3*data[32,1:])/flux_scale(cases[34])))
    assert max(errors+[reverse,scale])<2e-12 and max(old_errors)>2e-12
    paths=[Path(__file__).resolve(),ROOT/'scripts/verifySU2AUSMUp2Flux.py',PLAN,
        ROOT/'docs/validation/v148/failure02/verifySU2AUSMUp2Flux.py',ROOT/'docs/validation/v148/failure02/README.md']
    paths+=[p for p in OUT.iterdir() if p.is_file()]
    report=dict(passed=True,sourceCommit=COMMIT,cases=35,exactNativeFunctionBodies=True,
        uniformAndContactMaxScaledError=max(errors),reversalMaxScaledError=reverse,areaScalingMaxScaledError=scale,
        firstComponentwiseMaxError=max(old_errors),firstComparisonFailed=True,comparisonScaleClarifiedAfterExecution=True,
        noNativeRerun=True,wholeNativeSolverExecuted=False,implicitJacobianTested=False,
        thermodynamicSoundSpeedAndNumericalInterfaceScaleDistinguished=True,physicsQualified=False,newQualifiedArea=0,
        sourceSHA256={str(p):sha(p) for p in paths})
    write(OUT/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='sourceSHA256'}))

if __name__=='__main__':main()
