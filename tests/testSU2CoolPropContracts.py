"""检查错误热物性/旧守恒初态/不相容边界被拒绝，并实际核验EOS超域拒绝。"""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2ThermalRuntime import validate_coolprop_manifest
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=ROOT/'outputs/v121_validation/contracts';assert not out.exists();out.mkdir()
    case=ROOT/'outputs/v121_validation/thermal20';cfg=(case/'case.cfg').read_text();data=json.loads((case/'case_manifest.json').read_text())
    validate_coolprop_manifest(data,cfg)
    cases={
      'old_energy_restart':cfg.replace('RESTART_SOL= NO','RESTART_SOL= YES'),
      'unsupported_flux':cfg.replace('CONV_NUM_METHOD_FLOW= ROE','CONV_NUM_METHOD_FLOW= AUSM'),
      'ordinary_farfield':cfg+'\nMARKER_FAR= (farfield)\n',
      'unadapted_turbulence':cfg.replace('KIND_TURB_MODEL= NONE','KIND_TURB_MODEL= SST'),
      'mixed_wall_boundary':cfg+'\nMARKER_HEATFLUX= (wall,0)\n',
      'unreviewed_conductivity':cfg.replace('CONDUCTIVITY_MODEL= CONSTANT_PRANDTL','CONDUCTIVITY_MODEL= COOLPROP'),
      'unreviewed_prandtl':cfg.replace('PRANDTL_LAM= 0.72','PRANDTL_LAM= 0.71')}
    checks={}
    for name,text in cases.items():
        try:validate_coolprop_manifest(data,text)
        except AssertionError as e:checks[name]=dict(rejected=True,reason=str(e))
        else:raise AssertionError(name+' should be rejected')
    exe=ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-v1/bin/RAH_CoolProp_State.exe'
    src=out/'above_eos_range.csv';src.write_text('mode,a,b\nPT,5000,2100\n')
    output=out/'rejected_properties.csv';start=time.monotonic()
    p=subprocess.run([str(exe),str(src),str(output)],cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20,creationflags=subprocess.CREATE_NO_WINDOW)
    (out/'range_rejection.log').write_bytes(p.stdout)
    assert p.returncode!=0 and b'outside declared HEOS Air EOS range' in p.stdout
    checks['eos_2100K']=dict(rejected=True,exitCode=p.returncode,elapsedSeconds=time.monotonic()-start,executed=True,
        command=[str(exe),str(src),str(output)],executableSHA256=sha(exe))
    report=dict(checks=checks,validNativeConfigAccepted=True,engineRerun=False,physicsQualified=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/su2ThermalRuntime.py',case/'case.cfg',case/'case_manifest.json',src]})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
