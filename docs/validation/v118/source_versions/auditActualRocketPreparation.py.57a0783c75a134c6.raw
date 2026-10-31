"""独立核对真实原面/精确来流/全部守恒初场，成功后才创建一次性求解请求。"""
from pathlib import Path
import hashlib,json,re,uuid
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v118_validation';CASE=OUT/'actual01';ENGINE=CASE/'engine'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    target=OUT/'preflight_audit.json';assert not target.exists()
    job=read(OUT/'prepare01_process.json');assert job['status']=='terminal' and job['exitCode']==0 and all(job['sourceStableAfter'].values())
    orig=ROOT/'outputs/v56_validation/state_validity/state_validity_actual_rocket.mat'
    assert sha(orig)=='83c9f203d1f83ae1f3825da16809e521eeeb8640c03b0ce11839c2b9bf94a6af'
    s=loadmat(str(orig),squeeze_me=True,struct_as_record=False);fs=s['fs'];V=s['result'].V;F=s['result'].F.astype(int)-1
    data=loadmat(str(CASE/'mesh/hybrid_volume.mat'));old=loadmat(str(ROOT/'outputs/v113_validation/actual02/hybrid_volume.mat'))
    P=data['vertices'];W=data['wedges'].astype(int)-1;T=data['tetrahedra'].astype(int)-1;nv=len(V)
    assert np.array_equal(V,data['sourceVertices']) and np.array_equal(F+1,data['sourceFaces'])
    assert np.array_equal(P[:7*nv],old['vertices'][:7*nv]) and len(F)==19752 and len(W)==18*len(F)
    from sys import path
    path.insert(0,str(ROOT/'scripts'));from meshSU2HybridVolume import wedge_jacobian_min
    jac=wedge_jacobian_min(P[W]);assert min(jac)>0
    p=P[T];tet=np.einsum('ij,ij->i',p[:,1]-p[:,0],np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0]))/6;assert min(tet)>0
    manifest=read(ENGINE/'case_manifest.json');f=manifest['flight']
    for key in ['Mach','alpha','beta','temperature','pressure']:assert f[key]==getattr(fs,key)
    assert f['wallTemperature']==s['cfg'].flight.Twall==300 and f['iterations']==600 and f['namedCheckpoints']
    U=np.array(manifest['expectedVelocity']);assert np.max(abs(U-fs.flowVelocityVector))<1e-10
    initial=np.loadtxt(str(ENGINE/'initial_flow.csv'),delimiter=',',skiprows=1)
    assert np.array_equal(initial[:,0],np.arange(len(P))) and np.array_equal(initial[:,1:4],P)
    rho=initial[:,4];velocity=initial[:,5:8]/rho[:,None];kinetic=.5*np.sum(velocity**2,axis=1)
    e=initial[:,8]/rho-kinetic;temp=e*.4/287.05287;pressure=.4*rho*e;cp=1.4*287.05287/.4
    assert np.isfinite(initial).all() and min(rho)>0 and min(e)>0
    assert np.max(abs(pressure-fs.pressure))<1e-9
    assert np.max(abs(velocity[:nv]))==0 and np.max(abs(temp[:nv]-300))<1e-10
    delta=manifest['initialization']['thicknessM'];distance=np.linalg.norm(P[:19*nv].reshape(19,nv,3)-V[None,:,:],axis=2)
    b=np.ones(len(P));eta=np.minimum(1,distance.ravel()/delta);b[:19*nv]=eta*eta*(3-2*eta)
    H=cp*temp+kinetic;wanted=cp*300+(cp*fs.temperature+.5*np.dot(U,U)-cp*300)*b
    assert np.max(abs(H-wanted))/max(wanted)<1e-14
    assert np.max(abs(velocity-b[:,None]*U))<1e-10
    far=data['farfield'].astype(int).ravel()-1;assert np.all(b[far]==1) and np.all(b[18*nv:]==1)
    assert min(temp)>=min(300,fs.temperature)-1e-10
    assert np.max(abs(temp[18*nv:]-fs.temperature))<1e-10
    av=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2;area=float(np.linalg.norm(av,axis=1).sum())
    assert abs(area-1.9646292565715553)<1e-13
    config=(ENGINE/'case.cfg').read_text(encoding='utf-8')
    for line in ['CONV_RESIDUAL_MINVAL= -8','SOLVER= NAVIER_STOKES','KIND_TURB_MODEL= NONE','MGLEVEL= 0','SOLUTION_FILENAME= initial_flow']:
        assert line in config
    files=[Path(__file__).resolve(),orig,CASE/'mesh/hybrid_volume.mat',CASE/'mesh/mesh_process.json',CASE/'preparation.json',
           ENGINE/'case_manifest.json',ENGINE/'case.cfg',ENGINE/'mesh.su2',ENGINE/'projection.mat',ENGINE/'initial_flow.csv',
           ROOT/'writeSU2HybridInitialState.m',ROOT/'tests/prepareActualRocketFullDomain.m',ROOT/'scripts/meshSU2HybridVolume.py']
    report=dict(passed=True,originalFaces=len(F),areaM2=area,nodes=len(P),cells=len(W)+len(T),firstSixLayersExact=True,
        minimumJacobian=float(min(jac)),minimumTetVolume=float(min(tet)),temperatureRangeK=[float(min(temp)),float(max(temp))],
        totalEnthalpyRelativeError=float(np.max(abs(H-wanted))/max(wanted)),minimumDensity=float(min(rho)),
        exactOriginalFreestream=True,wallNoSlip=True,farfieldExact=True,initialGuessOnly=True,physicsQualified=False,
        inputSHA256={str(p):sha(p) for p in files})
    target.write_text(json.dumps(report,indent=2)+'\n')
    req=ENGINE/'engine_request.json';assert not req.exists()
    req.write_text(json.dumps(dict(workDirectory=str(ENGINE),configFile='case.cfg',meshFile='mesh.su2',threads=2,
        maximumSeconds=720,dryRun=False,jobToken=str(uuid.uuid4()),sourceFiles=[str(target),str(ENGINE/'projection.mat'),
        str(ENGINE/'case_manifest.json'),'writeSU2HybridInitialState.m','tests/prepareActualRocketFullDomain.m',
        'tests/auditActualRocketPreparation.py']),indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inputSHA256'},indent=2))
if __name__=='__main__':main()
