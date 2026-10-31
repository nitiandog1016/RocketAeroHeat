"""把同一完整体网格的定比热p/T/速度转换为HEOS初态，使用原生特征远场。"""
from pathlib import Path
import argparse,hashlib,json,os,re,sys
import numpy as np
import buildSU2ThermalRuntime as supervisor
ROOT=Path(__file__).resolve().parents[1]
RUNTIME=ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def setting(text,key):return re.search(r'^\s*'+key+r'\s*=\s*([^%\r\n]+)',text,re.M)[1].strip()
def prepare(folder,source,runtime_flavor='rah-coolprop-recon-v3',process_runner=None):
    assert runtime_flavor in ['rah-coolprop-v1','rah-coolprop-recon-v3']
    RUNTIME=ROOT/'tools/runtime/su2'/('8.6.0-'+runtime_flavor)
    folder=Path(folder).resolve();source=Path(source).resolve()
    for p in [folder,source]:p.relative_to(ROOT/'outputs')
    assert folder!=source and not (folder/'engine_process.json').exists()
    assert not (folder/'coolprop_initialization.json').exists() and not (folder/'initial_flow.csv').exists()
    cfg=(folder/'case.cfg').read_text();oldcfg=(source/'case.cfg').read_text()
    assert setting(cfg,'FLUID_MODEL')=='IDEAL_GAS' and setting(cfg,'RESTART_SOL')=='NO'
    assert setting(oldcfg,'FLUID_MODEL')=='IDEAL_GAS' and setting(oldcfg,'KIND_TURB_MODEL')=='NONE'
    r=read(source/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and not r['cancelled'] and all(r['inputsStable'].values())
    for name in ['mesh.su2','case.cfg']:assert sha(source/name)==r['inputSHA256'][str(source/name)]
    assert sha(folder/'mesh.su2')==sha(source/'mesh.su2'),'必须是同一原体网格'
    for key in ['MACH_NUMBER','AOA','SIDESLIP_ANGLE','FREESTREAM_TEMPERATURE','FREESTREAM_PRESSURE','MARKER_ISOTHERMAL','PRANDTL_LAM','VISCOSITY_MODEL','MU_REF','MU_T_REF','SUTHERLAND_CONSTANT']:
        assert setting(cfg,key)==setting(oldcfg,key),'原来流/热边界不同：'+key
    assert sha(source/'restart_flow.csv')==r['outputSHA256']['restart_flow.csv']
    manifest=read(folder/'case_manifest.json');assert manifest['thermalBoundary']['mode']=='isothermal'
    runtime=read(RUNTIME/'runtime.json');exe=Path(runtime['propertyExecutable'])
    assert sha(exe)==next(q['sha256'] for q in runtime['files'] if q['path']=='bin/RAH_CoolProp_State.exe')
    values=np.loadtxt(str(source/'restart_flow.csv'),delimiter=',',skiprows=1)
    assert values.shape[1]==9 and np.array_equal(values[:,0],np.arange(len(values))) and np.isfinite(values).all()
    rho=values[:,4];u=values[:,5:8]/rho[:,None];e=values[:,8]/rho-.5*np.sum(u*u,axis=1)
    R=float(setting(oldcfg,'GAS_CONSTANT'));gamma=float(setting(oldcfg,'GAMMA_VALUE'))
    T=e*(gamma-1)/R;p=rho*R*T;assert min(rho)>0 and min(T)>0 and min(p)>0
    out=folder/'initial_thermodynamics';assert not out.exists();out.mkdir()
    src=out/'PT_input.csv';dest=out/'PT_output.csv'
    pInf=float(setting(cfg,'FREESTREAM_PRESSURE'));TInf=float(setting(cfg,'FREESTREAM_TEMPERATURE'))
    with src.open('w') as f:
        f.write('mode,a,b\n');f.writelines('PT,%.17g,%.17g\n'%q for q in zip(p,T));f.write('PT,%.17g,%.17g\n'%(pInf,TInf))
    locked={str(q):sha(q) for q in [Path(__file__).resolve(),source/'case.cfg',source/'engine_process.json',source/'restart_flow.csv',folder/'mesh.su2',RUNTIME/'runtime.json',exe,src]}
    supervisor.OUT=out
    run=supervisor.run if process_runner is None else lambda *a:process_runner(out,*a)
    run('PT',[str(exe),str(src),str(dest)],ROOT,dict(os.environ,OMP_NUM_THREADS='1'),120)
    data=np.genfromtxt(str(dest),names=True,delimiter=',');assert np.array_equal(data['row'],np.arange(len(values)+1))
    local=data[:-1];freestream=data[-1];newrho=local['rho_kg_m3'];newe=local['e_J_kg']
    new=np.column_stack([values[:,:4],newrho,newrho[:,None]*u,newrho*(newe+np.sum(u*u,axis=1)/2)])
    initial=folder/'initial_flow.csv'
    np.savetxt(str(initial),new,delimiter=',',fmt=['%d']+['%.17e']*8,comments='',
        header='"PointID","x","y","z","Density","Momentum_x","Momentum_y","Momentum_z","Energy"')
    de=out/'DE_input.csv';deout=out/'DE_output.csv'
    # 从已写出的真实守恒CSV反算，验证写出舍入及动能扣除，而非仅查PT输出本身。
    written=np.loadtxt(str(initial),delimiter=',',skiprows=1);uu=written[:,5:8]/written[:,4,None]
    ee=written[:,8]/written[:,4]-.5*np.sum(uu*uu,axis=1)
    de.write_text('mode,a,b\n'+''.join('DE,%.17g,%.17g\n'%q for q in zip(written[:,4],ee)))
    run('DE',[str(exe),str(de),str(deout)],ROOT,dict(os.environ,OMP_NUM_THREADS='1'),120)
    check=np.genfromtxt(str(deout),names=True,delimiter=',')
    errors=dict(pressure=float(max(abs(check['p_Pa']/p-1))),temperature=float(max(abs(check['T_K']/T-1))),
        velocity=float(np.max(abs(uu-u)/np.maximum(abs(u),1))))
    assert max(errors.values())<1e-8 and all(sha(Path(q))==h for q,h in locked.items())
    speed=float(setting(cfg,'MACH_NUMBER'))*np.sqrt(gamma*R*TInf)
    alpha=np.deg2rad(float(setting(cfg,'AOA')));beta=np.deg2rad(float(setting(cfg,'SIDESLIP_ANGLE')))
    direction=np.array([np.cos(alpha)*np.cos(beta),np.sin(beta),np.sin(alpha)*np.cos(beta)])
    assert np.linalg.norm(direction)-1<1e-12
    mach=speed/freestream['a_m_s'];machvector=mach*direction
    for key,value in [('FLUID_MODEL','COOLPROP'),('CONV_NUM_METHOD_FLOW','ROE'),('MACH_NUMBER','%.17g'%mach),('RESTART_SOL','YES')]:
        cfg,n=re.subn(r'^\s*'+key+r'\s*=.*$',key+'= '+value,cfg,flags=re.M);assert n==1
    cfg,n=re.subn(r'^MARKER_FAR=.*\n','',cfg,flags=re.M);assert n==1
    cfg+='\nFLUID_NAME= Air\nMARKER_RIEMANN= (farfield, STATIC_SUPERSONIC_INFLOW_PT, %.17g, %.17g, %.17g, %.17g, %.17g)\n'%(pInf,TInf,*machvector)
    cfg+='SOLUTION_FILENAME= initial_flow\nREAD_BINARY_RESTART= NO\n'
    assert cfg.count('SOLUTION_FILENAME=')==1 and cfg.count('READ_BINARY_RESTART=')==1
    (folder/'case.cfg').write_text(cfg,encoding='utf-8')
    oldLineage=read(source/'continuation.json') if (source/'continuation.json').exists() else {}
    saved=[int(re.search(r'_(\d+)\.csv$',k)[1]) for k,h in r['outputSHA256'].items() if re.fullmatch(r'restart_flow_\d+\.csv',k) and h==r['outputSHA256']['restart_flow.csv']]
    assert saved
    lineage=dict(kind='same-mesh-primitive-to-HEOS-Air-initial-guess',thermodynamicModel='CoolProp-7.2.0-HEOS-Air',
        sourceDirectory=str(source),sourceLocalIteration=max(saved),sourceCumulativeIteration=oldLineage.get('completedUpdatesBeforeRun',0)+max(saved),
        samePhysicsAndNumerics=False,sameMesh=True,strictlyConservativeRemap=False,newPhysicsLocalIterationStartsAt=0,
        preserves=['p','T','velocity vector','mesh','wall boundary'],changes=['EOS density','EOS internal-energy reference','total conservative energy','convective flux','open boundary implementation'],
        nodeCount=len(values),errors=errors,sourceTemperatureRangeK=[float(min(T)),float(max(T))],
        originalVelocityMps=speed,originalReferenceMach=float(setting(oldcfg,'MACH_NUMBER')),EOSMach=float(mach),expectedVelocity=(speed*direction).tolist(),
        thermodynamicQuerySHA256={q.name:sha(q) for q in [src,dest,de,deout,out/'PT_process.json',out/'DE_process.json']},
        runtimeSHA256=sha(RUNTIME/'runtime.json'),initialStateSHA256=sha(initial),sourceSHA256=locked,sourceStable=True,
        maximumRelativeDensityChange=float(max(abs(newrho/rho-1))),maximumInternalEnergyDifferenceJkg=float(max(abs(newe-e))),physicsQualified=False)
    (folder/'coolprop_initialization.json').write_text(json.dumps(lineage,indent=2)+'\n')
    manifest['requiredRuntimeFlavor']=runtime_flavor;manifest['thermodynamicModel']=lineage['thermodynamicModel']
    manifest['thermalPropertyMode']='HEOS Air; Sutherland; k=mu*cp(T,rho)/Pr0.72'
    manifest['initialization']=lineage;manifest['expectedVelocity']=lineage['expectedVelocity']
    manifest['flight']['referenceMach']=manifest['flight']['Mach'];manifest['flight']['Mach']=float(mach)
    manifest['coolPropInitializationSHA256']=sha(folder/'coolprop_initialization.json')
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return lineage
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--folder',required=True);parser.add_argument('--source',required=True)
    parser.add_argument('--runtime-flavor',default='rah-coolprop-recon-v3',choices=['rah-coolprop-v1','rah-coolprop-recon-v3'])
    parser.add_argument('--verify-only',action='store_true');args=parser.parse_args()
    if args.verify_only:
        from su2ThermalRuntime import resolve
        folder=Path(args.folder).resolve();folder.relative_to(ROOT/'outputs')
        _,_,flavor,_,_=resolve(folder,(folder/'case.cfg').read_text())
        assert flavor==args.runtime_flavor
        r=read(folder/'coolprop_initialization.json');assert Path(r['sourceDirectory']).resolve()==Path(args.source).resolve()
    else:r=prepare(args.folder,args.source,args.runtime_flavor)
    print(json.dumps({k:r[k] for k in ['nodeCount','errors','sourceCumulativeIteration','originalVelocityMps','EOSMach','maximumRelativeDensityChange','maximumInternalEnergyDifferenceJkg']},indent=2))
