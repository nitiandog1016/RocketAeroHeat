"""真实早期原生场的来源/异常拒绝与能量零点核验，不运行新CFD。"""
from pathlib import Path
import copy,json,sys,time,shutil,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from readOpenFoamRocketState import BASE,read_json,read_frame
from openFoamAirThermo import HF


def main():
    start=time.monotonic();out=BASE/'state_tests01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_test.py')
    r=read_json(BASE/'rocket_run01/report.json');a=read_json(BASE/'rocket_run01/admission.json');name=next(iter(r['saved']));t=float(name)
    meta,f,g,m,ids,owners=read_frame(name,t,diagnostic=True)
    assert meta['status']=='interrupted-run-diagnostic' and meta['timeIndex']==1
    assert meta['wallHeatFluxWm2'] is None and meta['solidWallTemperatureK'] is None and not meta['engineeringHeatReady']
    assert np.array_equal(f['engineeringInternalEnergyJkg'],f['nativeSensibleInternalEnergyJkg']+HF)
    assert np.allclose(f['enthalpyJkg'],f['engineeringInternalEnergyJkg']+f['pressurePa']/f['densityKgm3'],rtol=0,atol=0)
    tri=g['points'][g['faces'][ids]];lengths=np.linalg.norm(tri-np.roll(tri,1,axis=1),axis=2);s=lengths.sum(axis=1)/2
    heron=np.sqrt(s*np.prod(s[:,None]-lengths,axis=1));area_error=float(abs(heron-g['area'][ids]).max())
    assert area_error<1e-14
    rows=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError) as e:rows.append(dict(name=name,rejected=True,reason=str(e)));return
        raise AssertionError('未拒绝 '+name)
    reject('failed_run_default',lambda:read_frame(name,t))
    reject('initial_field_as_solution',lambda:read_frame('0',0,True))
    reject('wrong_requested_time',lambda:read_frame(name,.00045,True))
    edits=[('wrong_profile',lambda r,a:r.update(profile='other')),
        ('missing_saved_frame',lambda r,a:r.update(saved={})),('wrong_source_time',lambda r,a:r['saved'][name].update(timeSeconds=.00045)),
        ('missing_energy',lambda r,a:r['saved'][name]['fieldsSHA256'].pop('e')),
        ('changed_field_sha',lambda r,a:r['saved'][name]['fieldsSHA256'].update(p='0'*64)),
        ('input_changed',lambda r,a:r['process']['inputStable'].update(anyInput=False)),
        ('live_group',lambda r,a:r['process'].update(liveGroupAfterCleanup={'123':{}})),
        ('log_unverified',lambda r,a:r['process'].update(logExportPassed=False)),
        ('admission_changed',lambda r,a:r.update(admissionSHA256='0'*64)),
        ('wrong_case',lambda r,a:a.update(case='outputs/other')),
        ('wrong_condition',lambda r,a:a['condition'].update(Mach=3)),
        ('wrong_CEA_speed',lambda r,a:a.update(velocityMps=[-1418.2630866,0,516.2055479]))]
    for label,edit in edits:
        rr=copy.deepcopy(r);aa=copy.deepcopy(a);edit(rr,aa)
        reject(label,lambda rr=rr,aa=aa:read_frame(name,t,True,rr,aa))
    p=subprocess.run([sys.executable,str(ROOT/'scripts/runOpenFoamRocketDevelopment.py')],capture_output=True,timeout=10)
    assert p.returncode and b'ValueError' in p.stderr and b'ModuleNotFoundError' not in p.stderr
    rows.append(dict(name='failed_candidate_cli_closed',rejected=True))
    report=dict(passed=True,positiveChecks=5,independentHeronAreaMaxErrorM2=area_error,rejections=rows,
        rejectionCount=len(rows),physicalQualification=False,CFDRuns=0,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['passed','positiveChecks','rejectionCount','elapsedSeconds']}))


if __name__=='__main__':main()
