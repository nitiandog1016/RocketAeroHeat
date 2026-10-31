"""核对已下载TMR原数据的形状与头部归一化；只确认参考输入，不宣布实验通过。"""
from pathlib import Path
import json,re,hashlib,io
import numpy as np
ROOT=Path(__file__).resolve().parents[1];DIR=ROOT/'docs/validation/v119/research'

def main():
    target=DIR/'reference_checks.json';assert not target.exists();reports=[]
    for name in ['M6Tw025','M8Tw048']:
        p=DIR/(name+'_Stat.dat');raw=p.read_bytes();meta=json.loads(p.with_suffix('.source.json').read_text())
        assert hashlib.sha256(raw).hexdigest()==meta['sha256']
        text=raw.decode('utf-8');rows=text.splitlines()
        i=next(i for i,s in enumerate(rows) if s.startswith('# Minf'))
        Mach,U,rho,T,Tw,ratio,delta_i=map(float,rows[i+1][1:].split())
        start=next(i for i,s in enumerate(rows) if s.startswith('#variables='))
        names=rows[start].split('=',1)[1].strip().split(',')
        a=np.loadtxt(io.StringIO('\n'.join(rows[start+1:])))
        count=int(re.search(r'(?m)^I=(\d+),',text)[1])
        assert a.shape==(count,len(names)) and a.shape[1]==28 and np.isfinite(a).all()
        assert np.all(np.diff(a[:,0])>0) and a[0,0]==0 and a[0,4]==0 and a[0,7]==1
        # 头部仅有限位数，以下是输入一致性容差，不是物理误差门槛。
        terr=abs(a[0,6]*T/Tw-1);merr=abs(U/(Mach*np.sqrt(1.4*287*T))-1)
        checks=dict(normalizedWallTemperatureConsistent=bool(terr<.003),headerMachConsistent=bool(merr<.003))
        reports.append(dict(case=name,Mach=Mach,speedMs=U,densityKgM3=rho,temperatureK=T,givenWallK=Tw,
                            Prandtl=.71,gasConstant=287,SutherlandCoefficient=1.458e-6,SutherlandTemperatureK=110.4,
                            profileRows=count,columns=len(names),normalizedWallTemperatureRoundingDifference=terr,
                            headerMachRoundingDifference=merr,inputChecks=checks,inputChecksPassed=all(checks.values()),
                            role='DNS numerical reference; no RocketAeroHeat flow comparison performed',
                            lowTemperatureModelMustMatch=True,experimentalValidation=False))
    target.write_text(json.dumps(dict(auditCompleted=True,allInputChecksPassed=all(r['inputChecksPassed'] for r in reports),
        toleranceUnchanged=.003,references=reports,sourceSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    print(json.dumps(reports,indent=2))

if __name__=='__main__':main()
