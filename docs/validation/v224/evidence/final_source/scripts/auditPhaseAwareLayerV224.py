"""模块职责：用Python独立读取MATLAB热链结果，核对新相态接口和解析平板，不重跑CFD/物性网格。"""
from pathlib import Path
import json,sys,subprocess
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=ROOT/'outputs/v224_validation';DOC=ROOT/'docs/validation/v224'

def main():
    out=BASE/'audit01';out.mkdir(exist_ok=False);source=BASE/'matlab01'
    assert read(BASE/'process_matlab01/report.json')['passed'] and read(source/'report.json')['passed']
    air=PhaseAwareViscousAir();records=[]
    # 独立语言、独立有界反算/梯形累积；这仍是关联式复算，不是实验验证。
    for name in ('uniform','variable','negative'):
        a=np.loadtxt(str(source/(name+'.csv')),delimiter=',',skiprows=1)
        s,p,te,u,tw=a[:,:5].T
        assert len(s)==81 and s[0]==0 and np.all(np.diff(s)>0)
        edge=air.evaluate(p,te);wall=air.evaluate(p,tw)
        tr=air.interpolated_transport(p,te)
        pr_edge=edge['cpJkgK']*tr['dynamicViscosityPas']/tr['conductivityWmK']
        haw=edge['enthalpyJkg']+np.sqrt(pr_edge)*u*u/2
        taw=air.inverse(p,haw,'h')
        href=.5*(edge['enthalpyJkg']+wall['enthalpyJkg'])+.22*(haw-edge['enthalpyJkg'])
        tref=air.inverse(p,href,'h');ref=air.evaluate(p,tref);transport=air.interpolated_transport(p,tref)
        pr=ref['cpJkgK']*transport['dynamicViscosityPas']/transport['conductivityWmK']
        nu=transport['dynamicViscosityPas']/ref['densityKgm3'];dh=haw-wall['enthalpyJkg']
        f=nu*u**5;integral=np.r_[0,np.cumsum(np.diff(s)*(f[:-1]+f[1:])/2)]
        theta=np.sqrt(.45*integral/u**6);cf=np.full(len(s),np.nan)
        for i in range(1,len(s)):
            j=max(0,i-4);lam=theta[i]**2/nu[i]*(u[i]-u[j])/(s[i]-s[j])
            assert -.09<lam<=.1
            cf[i]=2*(.2227+1.57*lam-1.8*lam**2)/(u[i]*theta[i]/nu[i])
        q=cf/(2*pr**(2/3))*edge['densityKgm3']*u*dh
        reference=np.column_stack([tref,href,taw,dh,pr,nu,theta,cf,q])
        # 原点热流/厚度未定义，不用解析theta=0覆盖NaN身份。
        error=np.max(abs(reference[1:]-a[1:,5:14])/np.maximum(1,abs(reference[1:])),axis=0)
        assert max(error)<3e-9,(name,error)
        if name=='uniform':
            blasius=.664/np.sqrt(u[1:]*s[1:]/nu[1:])
            exact_error=float(np.max(abs(cf[1:]/blasius-1)));assert exact_error<2e-4
        assert np.isnan(a[0,13]) and a[0,14]==8 and np.all(a[1:,14]==0)
        records.append(dict(case=name,maxScaledDifferences=error.tolist(),maxScaledDifference=float(max(error)),
            heatFluxRangeWm2=[float(min(q[1:])),float(max(q[1:]))],sourceSHA256=sha(source/(name+'.csv'))))
    # 现有前端只加载纯函数；SHA记录于回执，不复制并发代码。
    test=ROOT/'tests/testPhaseAwareLayerViewerV224.js'
    run=subprocess.run(['node',str(test)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15)
    (out/'frontend_test.txt').write_bytes(run.stdout);assert run.returncode==0,run.stdout
    save(out/'report.json',dict(passed=True,records=records,flatPlateBlasiusRelativeDifference=exact_error,
        propertyTableSHA256=sha(ROOT/'data/thermo/phase_aware_dry_air_v1.json'),
        matlabReportSHA256=sha(source/'report.json'),frontendSHA256=sha(ROOT/'web/viewer.js'),
        physicalQualification=False,newQualifiedAreaM2=0,browserInteractionTested=False))
    print(json.dumps(dict(passed=True,cases=3,maxScaledDifference=max(r['maxScaledDifference'] for r in records))))

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
        size=guard.size_bytes;guard.size_bytes=lambda unused:size(BASE)
        locks=[Path(__file__).resolve(),DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',
               ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
               ROOT/'data/thermo/phase_aware_dry_air_v1.json',ROOT/'web/viewer.js',
               ROOT/'tests/testPhaseAwareLayerViewerV224.js']
        locks += [BASE/'matlab01'/n for n in ('report.json','uniform.csv','variable.csv','negative.csv','workbench.json')]
        raise SystemExit(0 if guard.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],60,128,locks,ROOT,'offline') else 1)
