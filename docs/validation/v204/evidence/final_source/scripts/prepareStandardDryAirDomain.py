"""模块：复用已验证理想空气库生成新具名表；只读重用真实气体样本审核。"""
from pathlib import Path
import hashlib,json,time,shutil,sys
import numpy as np
from idealAirProperties import IdealAir
from nistIntegratedAir import AIR as OLD,R as OLD_R,X as OLD_X
from standardDryAir import StandardDryAir,MODEL,DOMAIN,pressure_ceiling
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def main():
    start=time.monotonic();out=ROOT/'outputs/v204_validation/standard01';out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    source=ROOT/'data/thermo/standard_dry_air_cp_v1.json';assert not source.exists()
    helper=IdealAir();legacy=json.loads((ROOT/'third_party/ideal_air/air_ideal.json').read_text(encoding='utf-8'))
    airref=json.loads((ROOT/'outputs/v204_validation/air02/report.json').read_text(encoding='utf-8'))
    assert airref['allPropertiesPassed'] and airref['allStableGasEvidence'] and airref['reference']=='Air'
    t=np.unique(np.r_[np.arange(100,2000.25,.5),298.15]);sample=helper.state('PT',100000.,t)
    baseline=helper.state('PT',100000.,298.15);hzero=float(baseline['h'][0])
    data=dict(schema='RAH.StandardDryAirCp.v1',model=MODEL,domain=DOMAIN,
        sourceModel=legacy['model'],sourceCommit=legacy['sourceCommit'],sourceSHA256=legacy['sourceSHA256'],
        sourceDOI=legacy['sourceDOI'],compositionSource='Lemmon et al.2000 Table3, p337',species=['N2','O2','Ar'],
        moleFractions=[.7812,.2096,.0092],gasConstantJkgK=legacy['gasConstant'],publishedMolarMassKgMol=legacy['molarMass'],
        temperatureK=t.tolist(),cpJkgK=sample['cp'].tolist(),oldOTHEnthalpyAt29815Jkg=hzero,
        reference='hs(298.15)=0; sf=0 relative entropy; hf=0',generatorSHA256=sha(__file__),
        legacyHelperBindings=helper.locks)
    save(source,data);model=StandardDryAir()
    # 模块一：新表连续积分/非节点采样与成熟理想项的独立关系核对。
    mid=(t[:-1]+t[1:])/2;query=np.r_[t[::19],mid[::17],100.,2000.]
    old=helper.state('PT',1000.,query);new=model.evaluate(1000.,query);scale=new['cpJkgK']*query
    tableErrors=dict(cp=float(np.max(abs(new['cpJkgK']/old['cp']-1))),
        h=float(np.max(abs(new['enthalpyJkg']-(old['h']-hzero))/scale)),
        e=float(np.max(abs(new['internalEnergyJkg']-(old['e']-hzero))/scale)),
        a=float(np.max(abs(new['soundSpeedMps']/old['sound']-1))))
    assert max(tableErrors.values())<=2e-6
    # 模块二：直接重用独立真实Air值，避免重跑相态或改变原误差门槛。
    aliases={'rho':'densityKgm3','a':'soundSpeedMps','cp':'cpJkgK','cv':'cvJkgK','h':'enthalpyJkg','e':'internalEnergyJkg'};rows=[]
    limits=airref['limits']
    for row in airref['rows']:
        p,T=row['p'],row['T'];d=model.evaluate(p,T);scale=float(d['cpJkgK']*T)
        errors={k:(float(d[key])-row['reference'][k])/(scale if k in {'h','e'} else row['reference'][k]) for k,key in aliases.items()}
        passed=all(abs(errors[k])<=limits[k] for k in limits)
        rows.append(dict(T=T,p=p,errors=errors,passed=passed,referencePhasePassed=row['stableGasEvidence']))
    maxima={k:max(abs(r['errors'][k]) for r in rows) for k in limits};passed=all(r['passed'] and r['referencePhasePassed'] for r in rows)
    assert passed
    # 模块三：原工况仅作输入设计；新旧模型差别显式报告，不制造发展场。
    pinf=1197.0279007383137;Tinf=226.50908361133003
    oldgas=OLD.evaluate(pinf,Tinf);newgas=model.evaluate(pinf,Tinf)
    differences={k:float(newgas[k]/oldgas[k]-1) for k in ['densityKgm3','soundSpeedMps','cpJkgK','cvJkgK']}
    names=['scripts/prepareStandardDryAirDomain.py','scripts/standardDryAir.py','docs/validation/v204/STANDARD_AIR_PLAN.md',
        'outputs/v204_validation/air02/report.json','data/thermo/standard_dry_air_cp_v1.json']
    (out/'mixture.dict').write_text(model.dictionary(),encoding='utf-8');save(out/'contract.json',model.contract())
    r=dict(passed=True,model=MODEL,domain=DOMAIN,tableNodes=len(t),tableTestStates=len(query),tableErrors=tableErrors,
        realEOSStates=len(rows),referenceMaxAbsoluteErrors=maxima,limits=limits,referenceSource='Lemmon2000 standardAir HEOS7.2.0',
        sourceScope='Paper standard-air model to2000K; high-temperature estimates partly based on nitrogen data, not uniformly experimental validation',
        realEOSGridPassed=passed,gridIsNotContinuousErrorProof=True,pythonDomainReady=True,nativeDomainDispatchVerified=False,
        actualCompositionChanged=True,oldMoleFractions=OLD_X.tolist(),newMoleFractions=data['moleFractions'],oldR=OLD_R,newR=model.R,
        freestreamRelativeChanges=differences,oldOTHEnthalpyAt29815Jkg=hzero,
        inputSHA256={n:sha(ROOT/n) for n in names},legacyHelperBindings=helper.locks,
        rows=rows,CFDRuns=0,builds=0,newQualifiedHeatAreaM2=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    assert time.monotonic()-start<60;save(out/'report.json',r)
    print(json.dumps({k:r[k] for k in ['passed','tableNodes','tableTestStates','tableErrors','realEOSStates','referenceMaxAbsoluteErrors','freestreamRelativeChanges','elapsedSeconds']}))

if __name__=='__main__':main()
