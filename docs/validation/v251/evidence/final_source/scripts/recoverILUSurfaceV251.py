"""模块职责：从更新前保存的合法旧六量恢复末次已完成状态的原面科研量，不补造原生末场。"""
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import brentq
import runNativeILUV251 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
from fasRecoveryFieldsV248 import binary
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC;CLOSED=True
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'

def recover():
    out=BASE/'recovery01';out.mkdir(exist_ok=False)
    audit=read(BASE/'audit01/report.json');refusal=audit['linearRefusal'];assert audit['passed'] and refusal
    step=refusal['oldStateIsCompletedUpdate'];assert step==audit['completedUpdates']==25
    assert not refusal['flowUpdateApplied'] and not refusal['SAUpdateApplied']
    a=binary(CASE/'linear_refused_state_direction.bin',b'RAH251RJ',335726,17);a=a[np.argsort(a[:,0])]
    assert np.array_equal(a[:,0],np.arange(len(a))) and np.isfinite(a).all()
    g={k:v for k,v in np.load(str(GEOMETRY)).items()};tri=g['triangles'];wall=g['wallIDs'];area=g['faceArea']
    rho=a[:,1];vel=a[:,2:5]/rho[:,None];e=a[:,5]/rho-.5*np.sum(vel**2,axis=1)-300000
    air=PhaseAwareViscousAir();T=inverse_energy(air,e);p=rho*air.R*T;phys=air.evaluate(p,T);transport=air.interpolated_transport(p,T)
    ids=np.unique(np.r_[np.linspace(0,len(a)-1,11).astype(int),np.argmin(T),np.argmax(T)])
    independent=np.array([brentq(lambda t:float(air.thermal(t)[1]-air.R*t)-e[i],60,2000,xtol=1e-10) for i in ids])
    terr=float(np.max(abs(T[ids]-independent)/np.maximum(1,abs(independent))));assert terr<2e-9
    # 法向差分是现有研究诊断定义；没有重算或冒称原生缓存壁面热流。
    q=np.full(len(a),np.nan);q[wall]=transport['conductivityWmK'][wall]*(T[g['normalNeighbor']]-300)/g['thermalDistance']
    xyz=g['points'][tri];lengths=[np.linalg.norm(xyz[:,i]-xyz[:,j],axis=1) for i,j in ((0,1),(1,2),(2,0))]
    semi=sum(lengths)/2;heron=np.sqrt(semi*(semi-lengths[0])*(semi-lengths[1])*(semi-lengths[2]))
    aerr=float(np.max(abs(heron-area)/area));assert aerr<2e-9 and abs(sum(area)-1.9646292565715553)<1e-12
    faces=np.column_stack((np.arange(1,19753),area,p[tri].mean(axis=1),q[tri].mean(axis=1),T[tri].mean(axis=1),np.full(19752,300)))
    assert np.isfinite(faces).all()
    target=out/('original_faces_%06d.csv'%step)
    np.savetxt(str(target),faces,delimiter=',',header='faceID,areaM2,pressurePa,heatFluxIntoWallWm2,wallFluidTemperatureK,givenWallTemperatureK',comments='',fmt='%.17g')
    restored=np.loadtxt(str(target),delimiter=',',skiprows=1);assert np.array_equal(restored,faces)
    report=dict(passed=True,update=step,sourceKind='pre-update-old-six-state-reconstruction',nativeFullFrame=False,nativeFinalCSV=False,
      terminalSpatialResidualAvailable=False,derivativeProbeUsed=False,states=len(a),BrentPoints=len(ids),maximumTemperatureRelativeError=terr,
      maximumIndependentHeronRelativeError=aerr,faces=19752,areaM2=float(sum(area)),faceCSV_SHA256=sha(target),sourceSHA256=sha(CASE/'linear_refused_state_direction.bin'),
      researchPowerIntoWallW=float(area@faces[:,3]),wallFluidMinimumK=float(min(T[wall])),wallFluidMaximumK=float(max(T[wall])),
      physicalQualification=False,solidTwSolved=False,newQualifiedAreaM2=0)
    save(out/'report.json',report);print({k:report[k] for k in ('passed','update','faces','researchPowerIntoWallW')})

def main():
    if CLOSED:raise SystemExit('v251 old-state surface recovery closed; follow NEXT_DECISION')
    locks=[Path(__file__).resolve(),Path(run.__file__),DOC/'SURFACE_RECOVERY.md',GEOMETRY,BASE/'audit01/report.json',CASE/'linear_refused_state_direction.bin',
      ROOT/'scripts/fasRecoveryFieldsV248.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'scripts/terminalCapacityV249.py']
    return supervise(run.setup_guard(),BASE,'recovery',[sys.executable,str(Path(__file__).resolve()),'task'],45,1536,locks,BASE,'audit')

if __name__=='__main__':
    if CLOSED:raise SystemExit('v251 old-state surface recovery closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':recover()
    else:raise SystemExit(0 if main() else 1)
