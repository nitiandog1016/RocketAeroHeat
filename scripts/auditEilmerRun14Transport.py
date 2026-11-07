"""独立比较Run14低密度输运；CoolProp仅为验证参考，不替换求解器或拟合热流。"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time


def audit(package,out):
    start=time.monotonic();out=Path(out);assert not out.exists();out.mkdir(parents=True)
    sys.path.insert(0,str(Path(package).resolve()))
    import CoolProp
    from CoolProp.CoolProp import PropsSI,get_global_param_string
    assert CoolProp.__version__=='7.2.0'
    # 官方N2数据库White(2006)固定常数；没有按实验q反标。
    R=8.31446261815324/.0280134;cp=3.5*R;rows=[]
    for T in [100.,120.4,150.,200.,295.2,500.,1000.,2000.]:
        mu=1.663e-5*(T/273.)**1.5*(273.+107.)/(T+107.)
        conductivity=mu*cp/.72
        refmu=PropsSI('V','T',T,'P',31.88,'Nitrogen')
        refk=PropsSI('L','T',T,'P',31.88,'Nitrogen')
        rows.append(dict(T=T,p=31.88,muSutherland=mu,kFixedPr=conductivity,
                         referenceMu=refmu,referenceK=refk,muRelativeDifference=mu/refmu-1,
                         kRelativeDifference=conductivity/refk-1))
    incoming=rows[1]
    # 事前筛选：100—200K黏度误差5%以内，密度误差不是输运误差。
    passed=all(abs(r['muRelativeDifference'])<=.05 for r in rows if r['T']<=200)
    rho=31.88/R/120.4;reference_rho=PropsSI('D','T',120.4,'P',31.88,'Nitrogen')
    h0=cp*120.4+2304**2/2
    ev=R*3393.44/math.expm1(3393.44/2467.)
    dlls={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(package).rglob('*.pyd')}
    report=dict(passed=passed,scope='Independent low-density transport correlation comparison; not experiment qualification',
        source='Lemmon & Jacobsen (2004), DOI 10.1023/B:IJOT.0000022327.04529.f3',
        CoolPropVersion=CoolProp.__version__,CoolPropGitRevision=get_global_param_string('gitrevision'),
        binarySHA256=dlls,rows=rows,rhoFromPT=rho,referenceRho=reference_rho,
        rhoVsReportedRelativeDifference=rho/.000881-1,gasConstant=R,cp=cp,
        perfectGasTotalEnthalpyJkg=h0,harmonicVibrationalEnergyJkg=ev,
        omittedVibrationFractionOfIncludingTotal=ev/(h0+ev),
        physicalQualification=False,temperatureAbove2000KReferenceChecked=False,
        elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8');return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--package',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();r=audit(a.package,a.out);print(json.dumps(r));raise SystemExit(0 if r['passed'] else 1)
