"""模块：只读绑定最终实现与已有数值；按已记录SHA恢复两个执行时源码快照。"""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
from standardDryAir import StandardDryAir,pressure_ceiling
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v204_validation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')

def main():
    out=BASE/'final_audit01';out.mkdir(exist_ok=False)
    current=(ROOT/'scripts/standardDryAir.py').read_text(encoding='utf-8')
    # 根据已记录变更恢复执行版，必须匹配当时原生监督SHA才写入；不是重新计算。
    before=current.replace("d=json.loads(SOURCE.read_text(encoding='utf-8'))\n        if d['model']!=MODEL or d['domain']!=DOMAIN:raise ValueError('STANDARD_AIR_MODEL_IDENTITY')", "d=json.loads(SOURCE.read_text(encoding='utf-8'));assert d['model']==MODEL and d['domain']==DOMAIN")
    before=before.replace("if (not np.isfinite(self.R) or self.R!=8.31451/.02896546 or self.T.shape!=self.Cp.shape\n                or self.T.ndim!=1 or len(self.T)<2 or not np.isfinite(self.T).all() or not np.isfinite(self.Cp).all()):\n            raise ValueError('STANDARD_AIR_INVALID_TABLE')\n        if self.T[0]!=100 or self.T[-1]!=2000 or not np.all(np.diff(self.T)>0) or not np.all(self.Cp>self.R):\n            raise ValueError('STANDARD_AIR_INVALID_TABLE')", "assert np.isfinite(self.R) and self.R>0 and self.T.shape==self.Cp.shape\n        assert self.T.ndim==1 and np.isfinite(self.T).all() and np.isfinite(self.Cp).all()\n        assert self.T[0]==100 and self.T[-1]==2000 and np.all(np.diff(self.T)>0) and np.all(self.Cp>self.R)")
    before=before.replace("i=np.flatnonzero(self.T==TREF)\n        if len(i)!=1:raise ValueError('STANDARD_AIR_MISSING_REFERENCE_NODE')\n        self.H-=self.H[i[0]];self.S-=self.S[i[0]]", "i=np.flatnonzero(self.T==TREF);assert len(i)==1;self.H-=self.H[i[0]];self.S-=self.S[i[0]]")
    admission=json.loads((BASE/'native01/admission.json').read_text(encoding='utf-8'));expected=admission['inputSHA256']['/mnt/e/RocketAeroHeat/scripts/standardDryAir.py']
    assert hashlib.sha256(before.encode('utf-8')).hexdigest()==expected
    (BASE/'native01/executed_standardDryAir.py').write_bytes(before.encode('utf-8'))
    test=(ROOT/'tests/testStandardDryAirDomain.py').read_text(encoding='utf-8').replace("out=ROOT/'outputs/v204_validation/domain_tests04'","out=ROOT/'outputs/v204_validation/domain_tests03'").replace("except ValueError:rejected.append(dict(case='wrong-'+key,message='identity rejected'))","except AssertionError:rejected.append(dict(case='wrong-'+key,message='identity rejected'))")
    report=json.loads((BASE/'domain_tests03/report.json').read_text(encoding='utf-8'));digest=report['sourceSHA256']['tests\\testStandardDryAirDomain.py']
    assert hashlib.sha256(test.encode('utf-8')).hexdigest()==digest;(BASE/'domain_tests03/executed_entry.py').write_bytes(test.encode('utf-8'))
    model=StandardDryAir();native=BASE/'native01';rows=np.array([[float(x) for x in s.split()[1:]] for s in (native/'positive/stdout.txt').read_text().splitlines() if s.startswith('ROW ')])
    d=model.evaluate(rows[:,0],rows[:,3]);keys=['cpJkgK','cvJkgK','enthalpyJkg','entropyJkgK','densityKgm3','soundSpeedMps']
    errors={k:float(np.max(abs(rows[:,i+5]-d[k])/np.maximum(1,abs(d[k])))) for i,k in enumerate(keys)}
    assert max(errors.values())<2e-9 and np.array_equal(rows[:,11],pressure_ceiling(rows[:,3]))
    for directory in ['library','probe']:
        original=ROOT/'native'/('openfoamStandardDryAir' if directory=='library' else 'openfoamStandardDryAirProbe')
        for path in original.rglob('*'):
            if path.is_file() and path.suffix in {'.C','.H'}:assert sha(path)==sha(native/directory/path.relative_to(original))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'));assert len(protected)==29
    assert all(sha(ROOT/p)==h for p,h in protected.items())
    # 实际API二元参数在交换CAS顺序后beta取倒数；核验不同模型的来源身份。
    source=json.loads((BASE/'source02/report.json').read_text(encoding='utf-8'))
    mix=json.loads((BASE/'mixture02/report.json').read_text(encoding='utf-8'))
    pairErrors=[]
    for actual in mix['binaryPairs']:
        row=next(x for x in source['binaryPairs'] if {x['CAS1'],x['CAS2']}=={actual['CAS1'],actual['CAS2']})
        for k in ['betaT','betaV','gammaT','gammaV','F']:
            value=row[k] if k not in ['betaT','betaV'] or row['CAS1']==actual['CAS1'] else 1/row[k]
            pairErrors.append(abs(float(actual[k])-value))
    assert max(pairErrors)<2e-15
    save(out/'report.json',dict(passed=True,protectedUnchanged=29,finalPythonNativeStates=17,propertyErrors=errors,
        originalPlanAcceptance=False,restoredExecutionSources=[dict(path='native01/executed_standardDryAir.py',sha256=expected),dict(path='domain_tests03/executed_entry.py',sha256=digest)],
        restoration='post-run source reconstruction verified against independent prior receipt SHA256',
        binaryPairMaxAbsDifference=max(pairErrors),finalPythonSHA256=sha(ROOT/'scripts/standardDryAir.py'),
        dataSHA256=sha(ROOT/'data/thermo/standard_dry_air_cp_v1.json'),CFDRuns=0,newQualifiedHeatAreaM2=0,physicalQualification=False))
    print(json.dumps(dict(passed=True,restoredSourceHashesMatched=True,protected=29,propertyErrors=errors)))

if __name__=='__main__':main()
