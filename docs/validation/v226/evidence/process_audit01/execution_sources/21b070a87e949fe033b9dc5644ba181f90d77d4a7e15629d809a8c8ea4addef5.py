"""模块职责：独立恢复原面外缘，核对指纹/面积/工作台并绘制已保存结果。"""
from pathlib import Path
import sys,json,hashlib,time,subprocess,re
import numpy as np
from scipy.optimize import brentq
from scipy.io import loadmat
from phaseAwareViscousAir import PhaseAwareViscousAir
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v226_validation'
def main():
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    source=BASE/'matlab02';assert read(BASE/'process_matlab02/report.json')['passed']
    result=read(source/'report.json');assert result['passed']
    e=read(source/'rocket_metadata.json');payload=read(source/'rocket_workbench.json')
    a=np.loadtxt(str(source/'rocket_faces.csv'),delimiter=',',skiprows=1);v=a[:,13].astype(bool);n=len(a)
    assert n==19752 and np.array_equal(a[:,0],np.arange(1,n+1))
    assert np.all(np.isnan(a[~v,2:11]))
    air=PhaseAwareViscousAir();cond=e['condition'];pinf=cond['pressurePa'];tinf=cond['temperatureK'];uinf=cond['speedMps']
    free=air.evaluate(pinf,tinf);h0=float(free['enthalpyJkg'])+.5*uinf*uinf
    T0=float(air.inverse(pinf,h0));shock=e['normalShock'];_,_,s2=air.thermal(shock['temperatureK']);_,_,s0=air.thermal(T0)
    p0=float(shock['pressurePa']*np.exp((s0-s2)/air.R));total=air.evaluate(p0,T0)
    theta=a[v,11];p=pinf+(p0-pinf)*np.cos(theta)**2
    wanted=s0+air.R*np.log(p/p0)
    te=np.array([brentq(lambda t:float(air.thermal(t)[2])-s,60,T0,xtol=1e-10) for s in wanted])
    state=air.evaluate(p,te);u=np.sqrt(2*(h0-state['enthalpyJkg']));transport=air.interpolated_transport(p,te)
    pre=state['cpJkgK']*transport['dynamicViscosityPas']/transport['conductivityWmK']
    taw=air.inverse(p,state['enthalpyJkg']+np.sqrt(pre)*u*u/2)
    reference=np.column_stack([p,te,u,state['densityKgm3'],u/state['soundSpeedMps'],taw])
    actual=a[v][:,[2,3,4,5,9,10]];difference=np.max(abs(reference-actual)/np.maximum(1,abs(reference)),axis=0)
    assert max(difference)<3e-9,difference
    # 原面几何/入射角及实际切向速度独立核对，旧压力与热值不读取。
    m=loadmat(str(ROOT/'outputs/v104_validation/baseline02/current_case.mat'),squeeze_me=True,struct_as_record=False)['current']
    V=m.V;F=m.F.astype(int)-1;c=V[F].mean(axis=1);cross=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])
    areas=np.linalg.norm(cross,axis=1)/2;normals=cross/(2*areas[:,None]);flow=np.array(cond['flowDirection'])
    radial=c-np.array(e['sphereCenter']);radial/=np.linalg.norm(radial,axis=1)[:,None]
    angles=np.arccos(np.clip(radial@(-flow),-1,1));assert max(abs(angles[v]-theta))<1e-12
    tangent=flow-(normals@flow)[:,None]*normals
    expectedU=u[:,None]*tangent[v]/np.linalg.norm(tangent[v],axis=1)[:,None]
    velocityError=float(np.max(abs(expectedU-a[v,6:9])));assert velocityError<1e-6
    assert np.max(abs(areas-a[:,1]))<1e-16 and abs(areas.sum()-1.9646292565715553)<2e-13
    assert np.array_equal(V,np.array(payload['vertices'])) and np.array_equal(F,np.array(payload['faces']))
    radius=e['geometryEstimate']['radius'];cap=max(0,1-e['geometryEstimate']['outerDepthRadiusFraction'])
    assert -flow[0]>=cap and np.all(radial[v,0]>=cap-1e-10) and np.all(theta<=np.deg2rad(cond['maximumPolarAngleDeg']))
    assert max(abs(a[v,12]-radius*theta))<1e-14
    grad=np.sqrt(2*(p0-pinf)/(float(total['densityKgm3'])*radius**2))
    assert abs(grad/e['stagnation']['velocityGradientPerSecond']-1)<1e-12
    # 与MATLAB相同字节合同独立打包；全部显示计算量的指纹不得来自自述accepted标记。
    digest=hashlib.sha256(e['sourceText'].encode('utf8'))
    arrays=[a[:,2],a[:,3],a[:,4],a[:,6:9],a[:,10],a[:,13],a[:,14],a[:,12],a[:,9],a[:,1]]
    for x in arrays:
        x=np.array(x,dtype='<f8');shape=x.shape if x.ndim==2 else (n,1)
        digest.update(np.array(shape,dtype='<f8').tobytes());digest.update(x.tobytes(order='F'))
    assert digest.hexdigest()==e['source']['sha256']
    meta=json.loads(e['sourceText']);assert meta['meshSHA256']==result['geometrySHA256']
    for name,value in meta['codeSHA256'].items():
        assert isinstance(value,str) and re.fullmatch('[0-9a-f]{64}',value)
        original=name[:-2]+'.m';assert sha(source/'executed_source'/original)==value
    assert meta['propertyTableSHA256']==sha(ROOT/'data/thermo/phase_aware_dry_air_v1.json')
    # 真实源JSON逐字段数值，未知热流仍是null，给定壁温仍是status2。
    for key,col in [('pressure',2),('edgeTemperature',3),('edgeSpeed',4),('edgeMach',9),('recoveryTemperature',10)]:
        f=next(f for f in payload['fields'] if f['key']==key);values=np.array(f['values'],float)
        assert np.array_equal(np.isfinite(values),v) and np.max(abs(values[v]-a[v,col]))<1e-9
    run=subprocess.run(['node',str(ROOT/'tests/testPhaseNoseViewerV226.js'),str(source)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend.txt').write_bytes(run.stdout);assert run.returncode==0,run.stdout
    plot(out,V,F,a,e)
    report=dict(passed=True,faces=n,candidateFaces=int(v.sum()),candidateAreaM2=float(areas[v].sum()),
      candidateAreaFraction=float(areas[v].sum()/areas.sum()),maxScaledDifferences=difference.tolist(),
      maximumScaledDifference=float(max(difference)),velocityMaximumAbsoluteErrorMps=velocityError,
      independentPayloadSHA256=digest.hexdigest(),sourceBindingVerified=True,
      finiteStagnationGradientPerSecond=float(grad),newHeatAreaM2=0,newQualifiedAreaM2=0,
      physicalQualification=False,browserInteractionTested=False,
      frontendSHA256=sha(ROOT/'web/viewer.js'),elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(report))
def plot(out,V,F,a,e):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import Normalize
    v=a[:,13].astype(bool);norm=Normalize(a[v,2].min()/1000,a[v,2].max()/1000);cmap=plt.get_cmap('viridis')
    colors=np.tile([.78,.79,.80,1.],(len(F),1));colors[v]=cmap(norm(a[v,2]/1000))
    fig=plt.figure(figsize=(14,5.3),facecolor='white')
    for j,zoom in enumerate((False,True)):
        ax=fig.add_subplot(1,3,j+1,projection='3d');mask=np.ones(len(F),bool)
        if zoom:mask=np.min(V[F,0],axis=1)>e['sphereCenter'][0]+.5*e['geometryEstimate']['radius']
        tri=V[F[mask]];coll=Poly3DCollection(tri,facecolors=colors[mask],edgecolors='none',rasterized=True);ax.add_collection3d(coll)
        pts=tri.reshape(-1,3);lo=pts.min(0);hi=pts.max(0);span=hi-lo
        for fun,a0,a1 in zip((ax.set_xlim,ax.set_ylim,ax.set_zlim),lo,hi):fun(a0,a1)
        ax.set_box_aspect(span);ax.view_init(elev=18,azim=-65);ax.set_xlabel('X [m]');ax.set_ylabel('Y [m]');ax.set_zlabel('Z [m]')
        ax.set_title('Original rocket: grey = unresolved' if not zoom else 'Nose pressure candidate [kPa]',fontsize=10)
        if zoom:
            bar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.035,pad=.10);bar.set_label('p [kPa]')
    ax=fig.add_subplot(1,3,3);angles=np.degrees(a[v,11]);order=np.argsort(angles)
    ax.plot(angles[order],a[v,2][order]/1000,color='#276c8e',lw=1.4,label='Pressure')
    ax.set_xlabel('Fitted-sphere polar angle [deg]');ax.set_ylabel('Pressure [kPa]',color='#276c8e');ax.grid(alpha=.2)
    right=ax.twinx();right.plot(angles[order],a[v,3][order],color='#a24e28',lw=1.4);right.set_ylabel('Edge temperature [K]',color='#a24e28')
    ax.set_title('Original face samples',fontsize=10)
    fig.suptitle('M5 / 30 km / alpha 20 deg | Phase-aware air | Engineering candidate; heat flux unresolved',fontsize=12)
    fig.subplots_adjust(left=.02,right=.94,top=.86,bottom=.17,wspace=.38)
    fig.text(.02,.025,'Original 19,752 faces; no interpolation into unresolved area. Pressure / temperature are not heat validation.',fontsize=9,color='#444444')
    fig.savefig(str(out/'rocket_nose_comparison.png'),dpi=160);plt.close(fig)
if __name__=='__main__':main()
