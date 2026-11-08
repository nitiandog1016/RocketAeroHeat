"""一次完整低温/相态/输运对照；复用旧热区证据，不扫描候选范围。"""
from pathlib import Path
import json,hashlib,time,shutil,sys,traceback
import numpy as np
from idealAirProperties import IdealAir
from phaseAwareDryAir import PhaseAwareDryAir,MODEL,DOMAIN,SOURCE,pressure_ceiling,dew_pressure
from standardDryAir import StandardDryAir
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v206_validation'
FIXED=ROOT/'tools/runtime/su2-thermal-build/source/subprojects/CoolProp'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))

def main():
    start=time.monotonic();out=BASE/'reference01';out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    rows=[]
    try:
        sys.path.insert(0,str(ROOT/'tools/runtime/coolprop-7.2.0'))
        import CoolProp as CP
        from CoolProp.CoolProp import AbstractState,get_global_param_string
        if CP.__version__!='7.2.0' or get_global_param_string('gitrevision')!='98b3523d5daa98454618d381d2ae53f7471d216b':raise ValueError('错误真实EOS版本')
        helper=IdealAir();old=StandardDryAir();original=read(FIXED/'dev/fluids/Air.json')
        # 模块一：旧表节点复用，仅生成此前缺少的60—100K节点。
        cold=np.arange(60.,100.,.5);cp=helper.state('PT',1.,cold)['cp'];data=dict(old.data)
        data.update(model=MODEL,domain=DOMAIN,temperatureK=np.r_[cold,old.T].tolist(),cpJkgK=np.r_[cp,old.Cp].tolist(),
            predecessorSHA256=sha(ROOT/'data/thermo/standard_dry_air_cp_v1.json'),
            phaseSource=original['ANCILLARIES']['pV'],transportSource=original['TRANSPORT'],
            fixedAirSourceSHA256=sha(FIXED/'dev/fluids/Air.json'),generatorSHA256=sha(Path(__file__).resolve()),legacyHelperBindings=helper.locks)
        if SOURCE.exists():raise ValueError('新表已存在；禁止自动覆写重试')
        save(SOURCE,data);gas=PhaseAwareDryAir()
        samples=np.r_[cold,cold+.25,100.];legacy=helper.state('PT',1.,samples);new=gas.evaluate(1.,samples)
        hzero=old.data['oldOTHEnthalpyAt29815Jkg'];scale=new['cpJkgK']*samples
        interpolation=dict(cp=float(np.max(abs(new['cpJkgK']/legacy['cp']-1))),
            h=float(np.max(abs(new['enthalpyJkg']-(legacy['h']-hzero))/scale)),
            e=float(np.max(abs(new['internalEnergyJkg']-(legacy['e']-hzero))/scale)),
            a=float(np.max(abs(new['soundSpeedMps']/legacy['sound']-1))))
        if max(interpolation.values())>2e-6:raise ValueError('新增表插值误差失败')
        # 模块二：只有新增低温区重新查询HEOS；旧318热区样本只读重用。
        state=AbstractState('HEOS','Air');state.update(CP.PT_INPUTS,1.,298.15);anchor=state.hmass()
        limits=dict(rho=.01,a=.01,cp=.02,cv=.02,h=.01,e=.01)
        aliases=dict(rho='densityKgm3',a='soundSpeedMps',cp='cpJkgK',cv='cvJkgK',h='enthalpyJkg',e='internalEnergyJkg')
        for t in np.arange(60.,100.,2.):
            for f in [.001,.01,.1,.5,1.]:
                if time.monotonic()-start>60:raise TimeoutError('物性总预算耗尽')
                p=float(pressure_ceiling(t))*f;state.unspecify_phase();state.update(CP.PT_INPUTS,p,t)
                ref=dict(rho=state.rhomass(),a=state.speed_sound(),cp=state.cpmass(),cv=state.cvmass(),h=state.hmass()-anchor,e=state.umass()-anchor)
                ideal=gas.evaluate(p,t);scale=float(ideal['cpJkgK']*t)
                err={k:(float(ideal[v])-ref[k])/(scale if k in ('h','e') else ref[k]) for k,v in aliases.items()}
                phase=int(state.phase());stable=phase in (int(CP.iphase_gas),int(CP.iphase_supercritical_gas)) and state.Q()<0 and state.isothermal_compressibility()>0
                transport=gas.transport(p,t);mu=state.viscosity();k=state.conductivity()
                row=dict(T=float(t),p=p,phase=phase,unforced=True,stable=bool(stable),reference=ref,errors=err,
                    viscosityRelativeError=float(transport['dynamicViscosityPas']/mu-1),conductivityRelativeError=float(transport['conductivityWmK']/k-1),
                    referenceViscosityPas=mu,referenceConductivityWmK=k,dewPressurePa=float(dew_pressure(t)))
                rows.append(row)
                if not stable or any(abs(err[k])>limits[k] for k in limits):raise ValueError('新增低温真实EOS/相态失败')
                if max(abs(row['viscosityRelativeError']),abs(row['conductivityRelativeError']))>.02:raise ValueError('低温输运差超过2%')
        prior=read(ROOT/'outputs/v204_validation/air02/report.json');replayed=[]
        for row in prior['rows']:
            current=gas.evaluate(row['p'],row['T']);scale=float(current['cpJkgK']*row['T'])
            errors={k:(float(current[v])-row['reference'][k])/(scale if k in ('h','e') else row['reference'][k]) for k,v in aliases.items()}
            if not row['stableGasEvidence'] or any(abs(errors[k])>limits[k] for k in limits):raise ValueError('旧318参照兼容失败')
            replayed.append(errors)
        # 模块三：新增输运范围对照；不是重复旧热物性网格。
        hot=[]
        for t in [100.,130.,180.,226.50908361133003,298.15,400.,600.,900.,1200.,1600.,2000.]:
            for factor in [.001,.1,1.]:
                p=float(pressure_ceiling(t))*factor;state.unspecify_phase();state.update(CP.PT_INPUTS,p,t);tr=gas.transport(p,t)
                row=dict(T=t,p=p,referenceViscosityPas=state.viscosity(),referenceConductivityWmK=state.conductivity(),
                    viscosityRelativeError=float(tr['dynamicViscosityPas']/state.viscosity()-1),conductivityRelativeError=float(tr['conductivityWmK']/state.conductivity()-1))
                hot.append(row)
                if max(abs(row['viscosityRelativeError']),abs(row['conductivityRelativeError']))>.02:raise ValueError('全温区输运差超过2%')
        names=['scripts/preparePhaseAwareDryAir.py','scripts/phaseAwareDryAir.py','data/thermo/phase_aware_dry_air_v1.json',
            'docs/validation/v206/PLAN.md','outputs/v204_validation/air02/report.json']
        files=[FIXED/'dev/fluids/Air.json',FIXED/'src/Backends/Helmholtz/TransportRoutines.cpp',FIXED/'LICENSE']
        for p in files:
            q=out/'fixed_source'/p.relative_to(FIXED);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        report=dict(passed=True,model=MODEL,domain=DOMAIN,oldModelUnchanged=True,newCpNodes=len(cold),totalCpNodes=len(gas.T),
            interpolation=interpolation,interpolationStates=len(samples),newEOSStates=len(rows),oldEOSRowsReplayed=len(replayed),newTransportHotStates=len(hot),
            coldMaxAbsoluteErrors={k:max(abs(r['errors'][k]) for r in rows) for k in limits},
            transportMaxRelativeErrors={k:max(abs(r[k]) for r in rows+hot) for k in ['viscosityRelativeError','conductivityRelativeError']},
            limits=limits,transportDifferenceLimit=.02,realEOSAnchorEnthalpyJkg=anchor,
            contract=gas.contract(),sourceSHA256={n:sha(ROOT/n) for n in names},fixedSourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in files},
            binarySHA256={p.name:sha(p) for p in (ROOT/'tools/runtime/coolprop-7.2.0/CoolProp').glob('*.pyd')},
            continuumQualified=False,finiteGridIsNotContinuousErrorProof=True,physicalQualification=False,CFDRuns=0,builds=0,
            coldRows=rows,hotTransportRows=hot,elapsedSeconds=time.monotonic()-start)
        save(out/'report.json',report);(out/'mixture.dict').write_text(gas.dictionary(),encoding='utf-8')
        print(json.dumps({k:report[k] for k in ['passed','newEOSStates','oldEOSRowsReplayed','coldMaxAbsoluteErrors','transportMaxRelativeErrors','elapsedSeconds']}))
    except BaseException as error:
        save(out/'report.json',dict(passed=False,error=repr(error),traceback=traceback.format_exc(),completedColdRows=rows,
            elapsedSeconds=time.monotonic()-start,CFDRuns=0,builds=0));raise

if __name__=='__main__':main()
