"""固定物性网格比较：相态失败和高温来源缺项显式保留，不改变生产压力门槛。"""
from pathlib import Path
import argparse,csv,hashlib,json,math,sys,time,shutil
import numpy as np
from nistIntegratedAir import AIR,R,X,PREF,TREF
ROOT=Path(__file__).resolve().parents[1]
LIMITS={'rho':.01,'a':.01,'cp':.02,'cv':.02,'h':.01,'e':.01}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def pressure_ceiling(t):return min(1e6,1e4*(float(t)/100)**3)

def ideal_state(p,t):
    """模块一：只用于审核的原理想关系重构，原生产evaluate压力合同原样保持。"""
    d={k:float(v) for k,v in AIR.evaluate(PREF,t).items()}
    return dict(rho=p/(R*t),a=d['soundSpeedMps'],cp=d['cpJkgK'],cv=d['cvJkgK'],h=d['enthalpyJkg'],e=d['internalEnergyJkg'])

def points():
    audit=json.loads((ROOT/'outputs/v203_validation/rocket_audit01/report.json').read_text(encoding='utf-8'))
    ts=np.unique(np.r_[np.arange(100,201,10),np.arange(200,1001,25),np.arange(1000,2001,100),
        298.15,499.999,500.001,699.999,700.001,999.999,1000.001,226.50908361133003,
        audit['conditionalTemperatureFromEnergyAtReferencePressureK'],audit['thermodynamicReference']['isentropicStagnationTemperatureK']])
    pts=[(float(t),pressure_ceiling(t)*f,'envelope') for t in ts for f in [.001,.01,.1,.5,1]]
    pts.extend([(226.50908361133003,1197.0279007383137,'freestream'),
        (audit['conditionalTemperatureFromEnergyAtReferencePressureK'],audit['rejectedInput']['pressurePa'],'v203-conditional-not-saved-field'),
        (audit['thermodynamicReference']['isentropicStagnationTemperatureK'],audit['thermodynamicReference']['isentropicStagnationPressurePa'],'isentropic-design-not-observed')])
    return pts

def run(kind,out):
    """模块二：原固定HEOS与实际相态求解；强制gas的诊断不参与稳定性验收。"""
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    sys.path.insert(0,str(ROOT/'tools/runtime/coolprop-7.2.0'))
    import CoolProp as CP
    from CoolProp.CoolProp import AbstractState,get_global_param_string,get_mixture_binary_pair_data
    assert CP.__version__=='7.2.0' and get_global_param_string('gitrevision')=='98b3523d5daa98454618d381d2ae53f7471d216b'
    model='Nitrogen&Oxygen&Argon' if kind=='mixture' else 'Air'
    state=AbstractState('HEOS',model)
    if kind=='mixture':state.set_mole_fractions(X.tolist())
    bindings={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'scripts/nistIntegratedAir.py',ROOT/'data/thermo/nist_air_shomate_v1.json',ROOT/'docs/validation/v204/REFERENCE_PLAN.md',
        ROOT/'outputs/v203_validation/rocket_audit01/report.json',*list((ROOT/'tools/runtime/coolprop-7.2.0/CoolProp').glob('*.pyd'))]}
    pairs=[]
    if kind=='mixture':
        for a,b in [('7727-37-9','7782-44-7'),('7727-37-9','7440-37-1'),('7782-44-7','7440-37-1')]:
            a,b=sorted((a,b))
            pairs.append(dict(CAS1=a,CAS2=b,**{k:get_mixture_binary_pair_data(a,b,k) for k in ['betaT','gammaT','betaV','gammaV','F','BibTeX']}))
    state.specify_phase(CP.iphase_gas);state.update(CP.PT_INPUTS,1.,TREF);hzero=state.hmass()
    anchor=dict(pressurePa=1,temperatureK=TREF,referenceEnthalpyJkg=hzero,Z=state.compressibility_factor(),forcedGas=True)
    rows=[]
    for t,p,label in points():
        assert time.monotonic()-start<60,'固定物性任务超过60秒'
        row=dict(T=t,p=p,label=label,ceiling=pressure_ceiling(t),kind=kind,unforcedPTPassed=False)
        state.unspecify_phase()
        try:
            state.update(CP.PT_INPUTS,p,t);row['unforcedPTPassed']=True
        except Exception as exc:
            row['unforcedError']=str(exc)
            try:state.specify_phase(CP.iphase_gas);state.update(CP.PT_INPUTS,p,t)
            except Exception as diagnostic:
                row.update(propertyPassed=False,diagnosticError=str(diagnostic));rows.append(row);continue
        q=state.Q();phase=int(state.phase());z=state.compressibility_factor();kt=state.isothermal_compressibility()
        ref=dict(rho=state.rhomass(),a=state.speed_sound(),cp=state.cpmass(),cv=state.cvmass(),h=state.hmass()-hzero,e=state.umass()-hzero)
        ideal=ideal_state(p,t);scale=ideal['cp']*t
        err={k:(ideal[k]-ref[k])/(scale if k in {'h','e'} else ref[k]) for k in ideal}
        passed=all(abs(err[k])<=LIMITS[k] for k in LIMITS)
        stable=bool(row['unforcedPTPassed'] and q<0 and kt>0 and ref['cp']>0 and ref['cv']>0 and .9<z<1.1)
        row.update(phaseCode=phase,Q=q,Z=z,isothermalCompressibilityPaInv=kt,reference=ref,ideal=ideal,errors=err,propertyPassed=passed,stableGasEvidence=stable)
        rows.append(row)
    # 模块三：实现一致、模型对照、来源适用和全物理资格分别登记。
    maxima={k:max(abs(r['errors'][k]) for r in rows if 'errors' in r) for k in LIMITS}
    report=dict(passed=True,kind=kind,reference=model,version=CP.__version__,revision=get_global_param_string('gitrevision'),
        compositionMoleFractions=state.get_mole_fractions(),gasConstantJmolK=state.gas_constant(),molarMassKgMol=state.molar_mass(),binaryPairs=pairs,
        anchor=anchor,limits=LIMITS,states=len(rows),propertyPassCount=sum(r['propertyPassed'] for r in rows),
        stableGasEvidenceCount=sum(r.get('stableGasEvidence',False) for r in rows),maxAbsoluteErrors=maxima,
        allPropertiesPassed=all(r['propertyPassed'] for r in rows),allStableGasEvidence=all(r.get('stableGasEvidence',False) for r in rows),
        exactComposition=(kind=='mixture'),wholeDomainExactMixtureSourceQualified=False,
        sourceCaveat='Gernert2013 actual mixture coefficients; no inherited Lemmon2000-mixture uncertainty or high-T guarantee' if kind=='mixture' else 'Lemmon2000 standard Air to2000K; different composition, not an exact-mixture qualification',
        productionAdmission=False,physicalQualification=False,CFDRuns=0,builds=0,inputSHA256=bindings,
        inputStable={p:sha(ROOT/p)==h for p,h in bindings.items()},rows=rows,elapsedSeconds=time.monotonic()-start)
    assert all(report['inputStable'].values());save(out/'report.json',report)
    with (out/'states.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['T_K','p_Pa','label','unforced_PT','stable_gas','property_pass',*['relative_'+k for k in LIMITS]])
        for row in rows:writer.writerow([row['T'],row['p'],row['label'],row['unforcedPTPassed'],row.get('stableGasEvidence',False),row['propertyPassed'],*[row.get('errors',{}).get(k,'') for k in LIMITS]])
    print(json.dumps({k:report[k] for k in ['kind','states','propertyPassCount','stableGasEvidenceCount','maxAbsoluteErrors','elapsedSeconds']}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--kind',choices=['mixture','air'],required=True);parser.add_argument('--out',required=True);args=parser.parse_args();run(args.kind,args.out)
