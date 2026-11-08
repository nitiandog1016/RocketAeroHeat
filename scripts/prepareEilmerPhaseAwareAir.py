"""模块职责：隔离模型登记与原生helper输入；不构建或运行流动任务。"""
from pathlib import Path
import json,hashlib,shutil,subprocess,time
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir
from phaseAwareDryAir import pressure_ceiling,SOURCE

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v210_validation'
COMMIT='f53f4609a0331d48efee69a4e4f3c3598378cc03'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def main():
    start=time.monotonic();out=BASE/'prepared01';out.mkdir(parents=True,exist_ok=False)
    src=ROOT/'outputs/v209_validation/source01/source/src/gas/init_gas_model.d'
    # 不修改共享factory，只增加一个精确的新模型分支。
    prior=json.loads((ROOT/'outputs/v209_validation/source01/report.json').read_text())
    expected=next(row['sha256'] for row in prior['sources'] if row['relativePath']=='src/gas/init_gas_model.d')
    assert sha(src)==expected and prior['sourceCommit']==COMMIT
    old=src.read_text(encoding='utf-8')
    assert old.count('import gas.gas_model;')==1 and old.count('case "IdealGas":')==1
    new=old.replace('import gas.gas_model;','import gas.gas_model;\nimport gas.rah_phase_aware_air;')
    new=new.replace('case "IdealGas":','case "rah-phase-aware-dry-air-eilmer-v1":\n            gm = new RAHPhaseAwareAir();\n            break;\n        case "IdealGas":')
    overlay=out/'overlay/gas';overlay.mkdir(parents=True)
    (overlay/'init_gas_model.d').write_text(new,encoding='utf-8')
    shutil.copyfile(ROOT/'native/eilmerPhaseAwareAir/gas/rah_phase_aware_air.d',overlay/'rah_phase_aware_air.d')
    shutil.copyfile(ROOT/'native/eilmerPhaseAwareAir/probe.d',out/'probe.d')
    (out/'gas.lua').write_text('model = "rah-phase-aware-dry-air-eilmer-v1"\n')
    (out/'invalid-gas.lua').write_text('model = "rah-phase-aware-dry-air-eilmer-v999"\n')
    model=PhaseAwareViscousAir()
    temps=[60.,60.25,80.,100.,132.5,180.,226.50908361133003,298.15,300.25,799.75,1000.,1274.25,1999.75,2000.]
    states=[];real=[];complex_rows=[];cases=[]
    def add(name,mode,x,y,dx=0.,dy=0.,expected='OK',meta=None):
        row='%s %s %.17g %.17g %.17g %.17g'%(name,mode,x,y,dx,dy)
        if dx==0 and dy==0:real.append(row)
        complex_rows.append(row)
        cases.append(dict(id=name,mode=mode,x=x if np.isfinite(x) else str(x),y=y if np.isfinite(y) else str(y),
            dx=dx,dy=dy,expected=expected,meta=meta))
    for i,t in enumerate(temps):
        p=float(min(1197.0279007383137,.2*pressure_ceiling(t)))
        st=model.evaluate(p,t);tr=model.interpolated_transport(p,t)
        vals=[p,t,float(st['densityKgm3']),float(st['internalEnergyJkg']),float(st['cpJkgK']),float(st['cvJkgK']),
              float(st['enthalpyJkg']),float(st['entropyJkgK']),float(st['soundSpeedMps']),float(tr['dynamicViscosityPas']),float(tr['conductivityWmK'])]
        states.append(dict(p=p,T=t,values=vals))
        for mode in ['PT','RT_RU','RT_RP','RT_RHOT','RT_PU','RT_PS','RT_HS']:
            add('state%d_%s'%(i,mode),mode,p,t,meta=dict(state=i))
        # 独立变量逐一扰动，不只验证恒等往返链。
        for mode,x,y in [('PT',p,t),('RU',vals[2],vals[3]),('RP',vals[2],p),('RHOT',vals[2],t),('PU',p,vals[3]),('PS',p,vals[7]),('HS',vals[6],vals[7])]:
            # 端点仅对原生同值正向/往返检查，不让跨实现端点一ulp差假称物理拒绝。
            if i in (0,len(temps)-1):continue
            for axis in (0,1):add('d%d_%s_%d'%(i,mode,axis),mode,x,y,1e-30 if axis==0 else 0.,1e-30 if axis==1 else 0.,meta=dict(state=i,derivativeAxis=axis))
    add('minimum','MIN',0,0)
    p60=float(pressure_ceiling(60.));p300=float(pressure_ceiling(300.))
    add('limit60','PT',p60,60.)
    add('limit300','PT',p300,300.)
    lo=float(model.thermal(60.)[1]-model.R*60.);hi=float(model.thermal(2000.)[1]-model.R*2000.)
    invalid=[('lowT','PT',1.,59.99),('highT','PT',1.,2000.01),('negativeP','PT',-1.,300.),('zeroP','PT',0.,300.),
        ('nanT','PT',1.,float('nan')),('infP','PT',float('inf'),300.),('dew','PT',p60*1.0001,60.),('pCeiling','PT',p300*1.0001,300.),
        ('lowE','RU',1e-5,lo-1),('highE','RU',1e-5,hi+1),('zeroRho','RU',0.,0.),('negativeRho','RP',-1.,1.),
        ('jointRU','RU',1.,lo+1),('jointRHOT','RHOT',1.,60.),('nanE','RU',1e-5,float('nan')),
        ('lowH','HS',float(model.thermal(60.)[1])-1,0.),('invalidPS','PS',-1.,0.),('badModel','BAD_MODEL',0.,0.)]
    for name,mode,x,y in invalid:add('bad_'+name,mode,x,y,expected='REJECT')
    for name,rows in [('real',real),('complex',complex_rows)]:(out/(name+'.tsv')).write_text('\n'.join(rows)+'\n')
    save(out/'cases.json',dict(states=states,cases=cases,columns=['p','T','rho','u','Cp','Cv','h','s','a','mu','k']))
    shutil.copyfile(__file__,out/'executed_prepare.py')
    save(out/'report.json',dict(passed=True,sourceCommit=COMMIT,factoryOriginalSHA256=sha(src),factoryOverlaySHA256=sha(overlay/'init_gas_model.d'),
        dataSHA256=sha(SOURCE),nativeSourceSHA256=sha(overlay/'rah_phase_aware_air.d'),baseStates=len(states),
        realRows=len(real),complexRows=len(complex_rows),elapsedSeconds=time.monotonic()-start,
        CFDRuns=0,physicalQualification=False))
    print(json.dumps(dict(states=len(states),realRows=len(real),complexRows=len(complex_rows),seconds=time.monotonic()-start)))

if __name__=='__main__':raise ValueError('v210气体准备已完成，保留准确历史输入；不重复准备。')
