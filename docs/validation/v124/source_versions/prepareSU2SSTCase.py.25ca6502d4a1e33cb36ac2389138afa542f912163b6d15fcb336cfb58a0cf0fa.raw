"""显式SST-2003m定比热研究入口；新模型初始化，不冒充同物理续算。"""
from pathlib import Path
import argparse,json,re
import numpy as np
from prepareSU2CheckpointContinuation import config,CONTROL,sha,read,ROOT
from su2SSTInterface import state,wall_geometry


def set_value(text,key,value):
    text,n=re.subn(r'^\s*'+re.escape(key)+r'\s*=.*$',key+'= '+str(value),text,flags=re.M)
    assert n<=1
    return text if n else text+'\n'+key+'= '+str(value)+'\n'


def prepare(folder,source=None,intensity=.01,viscosity_ratio=10.,prt=.9,provenance=''):
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    assert not (folder/'engine_process.json').exists() and not (folder/'sst_contract.json').exists()
    assert provenance.strip(), '必须说明全湍流和入口参数来源，不默认当作测量值'
    assert all(np.isfinite([intensity,viscosity_ratio,prt])) and 0<intensity<1 and viscosity_ratio>0 and prt>0
    b=config(folder/'case.cfg');manifest=read(folder/'case_manifest.json')
    assert manifest.get('requiredRuntimeFlavor','official')=='official'
    for k,v in {'SOLVER':'NAVIER_STOKES','KIND_TURB_MODEL':'NONE','FLUID_MODEL':'IDEAL_GAS',
                'SYSTEM_MEASUREMENTS':'SI','REF_DIMENSIONALIZATION':'DIMENSIONAL','MGLEVEL':'0',
                'VISCOSITY_MODEL':'SUTHERLAND','CONDUCTIVITY_MODEL':'CONSTANT_PRANDTL'}.items():assert b[k]==v,k
    assert b['RESTART_SOL']=='NO' and 'MARKER_ISOTHERMAL' in b and 'MARKER_WALL_FUNCTIONS' not in b
    assert b.get('TIME_DOMAIN','NO')=='NO'
    gamma=float(b['GAMMA_VALUE']);R=float(b['GAS_CONSTANT']);T=float(b['FREESTREAM_TEMPERATURE'])
    rho=float(b['FREESTREAM_PRESSURE'])/(R*T);speed=float(b['MACH_NUMBER'])*np.sqrt(gamma*R*T)
    def mu(t):return float(b['MU_REF'])*(t/float(b['MU_T_REF']))**1.5*(float(b['MU_T_REF'])+float(b['SUTHERLAND_CONSTANT']))/(t+float(b['SUTHERLAND_CONSTANT']))
    kinf=1.5*(speed*intensity)**2;omega=rho*kinf/(mu(T)*viscosity_ratio)
    contract=dict(schemaVersion=1,model='SST-2003m',assumption='fully-turbulent; transition unknown',
        intensity=intensity,viscosityRatio=viscosity_ratio,Prt=prt,provenance=provenance,
        kInfinity=kinf,omegaInfinity=omega,wallTreatment='resolved smooth wall; no wall function',
        thermodynamics='calorically-perfect ideal gas',turbulentConvection='first-order scalar upwind',
        compressibilityCorrection='none',physicsQualified=False,newQualifiedArea=0)
    text=(folder/'case.cfg').read_text(encoding='utf-8-sig')
    changes={'SOLVER':'RANS','KIND_TURB_MODEL':'SST','SST_OPTIONS':'V2003m',
        'FREESTREAM_TURBULENCEINTENSITY':format(intensity,'.17g'),'FREESTREAM_TURB2LAMVISCRATIO':format(viscosity_ratio,'.17g'),
        'TURBULENT_CONDUCTIVITY_MODEL':'CONSTANT_PRANDTL_TURB','PRANDTL_TURB':format(prt,'.17g'),
        'CONV_NUM_METHOD_TURB':'SCALAR_UPWIND','MUSCL_TURB':'NO','TIME_DISCRE_TURB':'EULER_IMPLICIT','CFL_REDUCTION_TURB':'1.0',
        'CONV_FIELD':'(RMS_DENSITY, RMS_ENERGY, RMS_TKE, RMS_DISSIPATION)',
        'WRT_RESTART_OVERWRITE':'NO','WRT_SURFACE_OVERWRITE':'NO'}
    if source is not None:
        source=Path(source).resolve();source.relative_to(ROOT/'outputs');assert source!=folder
        r=read(source/'engine_process.json');a=config(source/'case.cfg')
        assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and not r['cancelled']
        assert all(r['inputsStable'].values()) and r['runtimeFlavor']=='official'
        assert {k:v for k,v in a.items() if k not in CONTROL}=={k:v for k,v in b.items() if k not in CONTROL}, '源场与目标的原物理/数值输入不一致'
        for name in ['case.cfg','mesh.su2']:assert sha(source/name)==r['inputSHA256'][str(source/name)]
        assert sha(source/'mesh.su2')==sha(folder/'mesh.su2')
        for name in ['restart_flow.csv','surface_flow.vtk','history.csv']:assert sha(source/name)==r['outputSHA256'][name]
        assert a['WRT_RESTART_OVERWRITE']=='NO' and a['WRT_SURFACE_OVERWRITE']=='NO'
        pairs=[]
        for p in source.glob('restart_flow_*.csv'):
            match=re.fullmatch(r'restart_flow_(\d+)\.csv',p.name)
            if not match:continue
            w=source/('surface_flow_'+match[1]+'.vtk')
            if w.is_file() and sha(p)==r['outputSHA256']['restart_flow.csv'] and sha(w)==r['outputSHA256']['surface_flow.vtk']:
                assert sha(p)==r['outputSHA256'][p.name] and sha(w)==r['outputSHA256'][w.name];pairs.append(int(match[1]))
        assert pairs;step=max(pairs)
        hist=np.atleast_1d(np.genfromtxt(str(source/'history.csv'),names=True,delimiter=','))
        assert np.array_equal(hist['Inner_Iter'],np.arange(len(hist))) and step==int(hist['Inner_Iter'][-1])
        data,keys,nd,density,u,e=state(source/'restart_flow.csv',False)
        g=wall_geometry(folder/'mesh.su2');assert np.allclose(g['points'],data[:,1:nd+1],rtol=0,atol=1e-14)
        k=np.full(len(data),kinf);w=np.full(len(data),omega);ids=g['wallIDs'];k[ids]=0
        w[ids]=60*mu(e[ids]*(gamma-1)/R)/(density[ids]*.075*g['nearestWallDistance']**2)
        output=np.column_stack([data,k,w]);output[:,keys.index('Energy')]+=density*k
        path=folder/'initial_flow.csv';assert not path.exists()
        np.savetxt(str(path),output,delimiter=',',fmt='%.17g',header=','.join('"'+s+'"' for s in keys+['Turb_Kin_Energy','Omega']),comments='')
        _,_,_,_,_,back=state(path,True);error=float(max(abs(back-e)/e));assert error<1e-12
        locks=[source/n for n in ['engine_process.json','case.cfg','mesh.su2','restart_flow.csv','surface_flow.vtk','history.csv']]
        lineage=dict(kind='same-mesh-laminar-to-SST-initial-guess',sourceDirectory=str(source),sourceLocalIteration=step,
            sourceCumulativeIteration=None,newPhysicsLocalIterationStartsAt=0,sameMesh=True,samePhysicsAndNumerics=False,
            initialStateSHA256=sha(path),thermalEnergyRelativeError=error,initialGuessOnly=True,physicsQualified=False)
        if (source/'continuation.json').exists():
            assert sha(source/'continuation.json')==r['inputSHA256'][str(source/'continuation.json')]
            old=read(source/'continuation.json');lineage['sourceCumulativeIteration']=int(old['completedUpdatesBeforeRun'])+step
            locks.append(source/'continuation.json')
        if (source/'initial_flow.transfer.json').exists():
            assert sha(source/'initial_flow.transfer.json')==r['inputSHA256'][str(source/'initial_flow.transfer.json')]
            old=read(source/'initial_flow.transfer.json');lineage['upstreamBeforeMeshTransferIteration']=old['sourceCumulativeIteration']
            locks.append(source/'initial_flow.transfer.json')
        lineage['sourceSHA256']={str(p):sha(p) for p in locks}
        lineage['nearestWallDistanceRangeM']=[float(min(g['nearestWallDistance'])),float(max(g['nearestWallDistance']))]
        (folder/'sst_initialization.json').write_text(json.dumps(lineage,indent=2)+'\n',encoding='utf-8')
        contract['initializationSHA256']=sha(folder/'sst_initialization.json');manifest['initialization']=lineage
        changes.update(RESTART_SOL='YES',SOLUTION_FILENAME='initial_flow',READ_BINARY_RESTART='NO')
    else:contract['initialization']='native uniform freestream; wall BC imposed by solver'
    for k,v in changes.items():text=set_value(text,k,v)
    (folder/'case.cfg').write_text(text,encoding='utf-8')
    contract['configurationSHA256']=sha(folder/'case.cfg');contract['meshSHA256']=sha(folder/'mesh.su2')
    (folder/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    manifest['sst']=contract;manifest['sstContractSHA256']=sha(folder/'sst_contract.json')
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    return contract


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--source');p.add_argument('--intensity',type=float,default=.01)
    p.add_argument('--viscosity-ratio',type=float,default=10.);p.add_argument('--prt',type=float,default=.9);p.add_argument('--provenance',required=True);a=p.parse_args()
    r=prepare(a.folder,a.source,a.intensity,a.viscosity_ratio,a.prt,a.provenance);print(json.dumps(r,ensure_ascii=False))
