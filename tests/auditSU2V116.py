"""独立核对600步配对、原面导出、成组原图读数；不重算流场或赋予物理资格。"""
from pathlib import Path
import hashlib,json,csv,re,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v116_validation'
sys.path.insert(0,str(ROOT/'tests'))
from auditSU2V110 import vtk_tokens
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    target=OUT/'independent_audit.json';assert not target.exists()
    records={};paths=[Path(__file__).resolve(),ROOT/'tests/auditSU2V110.py'];surfaces={};histories={}
    for name in ['control','candidate']:
        folder=OUT/'study01'/name;ui=OUT/'review02';receipt=read(folder/'engine_process.json')
        assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut']
        assert all(receipt['inputsStable'].values()) and receipt['request']['threads']==2
        for f,h in receipt['inputSHA256'].items():assert sha(Path(f))==h,f
        for f,h in receipt['outputSHA256'].items():assert sha(folder/f)==h,f
        h=np.atleast_1d(np.genfromtxt(str(folder/'history.csv'),delimiter=',',names=True));histories[name]=h
        assert np.array_equal(h['Inner_Iter'],np.arange(600)) and (h['HF']==0).all()
        r=loadmat(str(ui/(name+'_result.mat')),squeeze_me=True,struct_as_record=False)['r']
        V=r.V;F=r.F.astype(int)-1;x,fields=vtk_tokens(folder/'surface_flow.vtk')
        distance,ids=cKDTree(V).query(x);assert max(distance)<1e-6 and len(np.unique(ids))==len(V)
        values=np.zeros((len(V),3));values[ids]=np.column_stack([fields[k][:,0] for k in ('Pressure','Heat_Flux','Temperature')])
        expected=values[F].mean(axis=1)
        for j,k in enumerate(['pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw']):assert np.allclose(expected[:,j],getattr(r,k),rtol=1e-12,atol=1e-8)
        assert not r.productionQualified and r.heatFluxIsBoundaryInput and np.isnan(r.wallTemperatureInput).all()
        assert r.convergenceAudit.terminalIteration==599 and not r.convergenceAudit.checkpointIdentityConfirmed
        table=np.genfromtxt(str(ui/(name+'_fields.csv')),delimiter=',',names=True)
        assert np.array_equal(table['originalFaceID'],np.arange(1,len(F)+1))
        for j,k in enumerate(['pressureRaw_Pa','heatFluxRaw_W_m2','fluidNodeTemperatureMean_K']):assert np.allclose(table[k],expected[:,j],rtol=1e-12,atol=1e-8)
        assert table['heatFluxIsBoundaryInput'].all() and not table['productionQualified'].any()
        assert np.isnan(table['givenWallTemperature_K']).all() and np.isnan(table['solvedWallTemperature_K']).all()
        with (ui/(name+'_section.csv')).open(encoding='utf-8-sig',newline='') as stream:rows=list(csv.reader(stream))[1:]
        for row in rows:
            assert abs(float(row[0])-.158496)<1e-14
            assert abs(float(row[10])-expected[int(row[2])-1,0])<1e-8
        assert len(rows)==32
        rest=np.loadtxt(str(folder/'restart_flow.csv'),delimiter=',',skiprows=1)
        rho=rest[:,4];kinetic=.5*np.sum((rest[:,5:8]/rho[:,None])**2,axis=1);e=rest[:,8]/rho-kinetic
        T=e*.4/287.05287;p=.4*rho*e;H=1.4*e+kinetic;H0=1.4*287.05287/.4*(710*5/9)
        assert np.isfinite(rest).all() and min(rho)>0 and min(e)>0
        # 同态EOS只能确认两个保存文件一致；HF全零不能推定其迭代身份。
        distance,node=cKDTree(rest[:,1:4]).query(V);assert max(distance)<1e-14
        eos=max(float(np.max(abs(p[node]/values[:,0]-1))),float(np.max(abs(T[node]/values[:,2]-1))))
        assert eos<2e-6
        text=(folder/'engine_stdout.log').read_text()
        assert 'All volume elements are correctly oriented.' in text and 'All surface elements are correctly oriented.' in text
        rec=read(OUT/'initial_comparison.json')[name]
        assert rec['nonphysicalReconstructionCount']==int(re.findall(r'Warning: (\d+) reconstructed states',text)[-1])
        assert np.allclose(rec['temperatureRange'],[min(T),max(T)],rtol=1e-12)
        rec.update(originalFaces=len(F),CSVRows=len(table),sectionRows=len(rows),maximumWallEOSError=eos,
            rhoResidualRatio=float(10**(h['rmsRho'][-1]-h['rmsRho'][0])),
            energyResidualRatio=float(10**(h['rmsRhoE'][-1]-h['rmsRhoE'][0])),
            savedIterationConfirmed=False,savedSurfaceAndRestartSameState=True,physicsQualified=False,newQualifiedArea=0)
        records[name]=rec;surfaces[name]=(V,F,values,expected)
        paths += [folder/n for n in ['case.cfg','engine_process.json','surface_flow.vtk','restart_flow.csv','history.csv']]
        paths += [ui/(name+s) for s in ['_fields.csv','_section.csv','_result.mat']]
    V,F,_,a=surfaces['control'];VV,FF,_,b=surfaces['candidate']
    assert np.array_equal(V,VV) and np.array_equal(F,FF)
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    difference=dict(pressureAreaWeightedRelativeL1=float(area@abs(a[:,0]-b[:,0])/(area@abs(a[:,0]))),
                    wallFluidTemperatureAreaWeightedRelativeL1=float(area@abs(a[:,2]-b[:,2])/(area@abs(a[:,2]))),
                    meaning='unconverged candidate difference, not discretization error or validation')
    # 用独立线性最小二乘逆式复算图11，不调用MATLAB读图程序。
    rawfile=ROOT/'tests/data/schaefer_ferguson_1961/heat_fig11.json';raw=read(rawfile);c=raw['calibration']
    source=ROOT/raw['archivedSourceImage'];assert sha(source)==raw['sourceImageSHA256']
    fit=np.linalg.lstsq(np.column_stack([np.log10(c['heatStSqrtRe']),np.ones(3)]),c['heatY'],rcond=None)[0]
    with (OUT/'review02/paired_heat_data.csv').open(encoding='utf-8-sig',newline='') as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==len(raw['points'])==40;matches=0
    for point,row in zip(raw['points'],rows):
        assert row['id']==point['id']
        if point['heatPixel'] is None:
            assert row['missing']=='1' and np.isnan(float(row['St']));continue
        expectedSt=10**((point['heatPixel'][1]-fit[1])/fit[0])/np.sqrt(240000)
        assert abs(float(row['St'])/expectedSt-1)<1e-12
        for key,pixel in [('Tw_K','wallPixel'),('Taw_K','adiabaticPixel')]:
            expectedT=(478-point[pixel][1])/277*(710*5/9)
            assert abs(float(row[key])-expectedT)<1e-10
        expectedUpper=172+553*point['station_in']/18;expectedLower=170+555*point['station_in']/18
        passed=max(abs(point['wallPixel'][0]-expectedUpper),abs(point['adiabaticPixel'][0]-expectedUpper),abs(point['heatPixel'][0]-expectedLower))<=1.5
        assert passed==(row['nominalStationMatched']=='1');matches+=int(passed)
    # 已保存2999场的近壁外端诊断：用真实节点列，不冒称边界层厚度。
    oldpath=ROOT/'outputs/embedded_workbench/96284313-167b-4afc-a198-39d8a5cc2385/restart_flow_002999.csv'
    old=np.loadtxt(str(oldpath),delimiter=',',skiprows=1);nv=len(V);profiles=[]
    for xx in [.33825,.291358,.236650,.189758,.142646]:
        i=int(np.argmin(abs(V[:,0]-xx)+abs(V[:,2])+abs(V[:,1]-.02286)))
        nodes=np.arange(25)*nv+i;d=old[nodes];vel=d[:,5:8]/d[:,4,None]
        temperature=(d[:,8]/d[:,4]-.5*np.sum(vel**2,axis=1))*.4/287.05287
        profiles.append(dict(wallNodeID=i,wallX=float(V[i,0]),distanceM=np.linalg.norm(d[:,1:4]-d[0,1:4],axis=1).tolist(),
            axialVelocityOverFreestream=(-vel[:,0]/812.1592303453607).tolist(),temperatureK=temperature.tolist()))
    (OUT/'initial_wall_profiles.json').write_text(json.dumps(dict(source=str(oldpath),sha256=sha(oldpath),profiles=profiles,convergedBoundaryLayerThickness=False),indent=2)+'\n')
    paths += [rawfile,source,oldpath,OUT/'review02/paired_heat_data.csv',OUT/'initial_comparison.json']
    report=dict(cases=records,candidateDifference=difference,pairedRecords=40,missingRecords=4,
        brokenRecords=4,nominalStationMatchedRecords=matches,digitizationStatus='draft; no point accepted for experimental comparison',
        candidateAdopted=False,physicalSolutionQualified=False,newQualifiedArea=0,
        inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths})
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    # 科研静态图：只比较已保存研究场，明确绝热给定和未收敛，避免绘制虚假的热流验证。
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10.5,7.0))
    for name,color in [('control','#345a7d'),('candidate','#a46533')]:
        h=histories[name];v,f,nodal,_=surfaces[name]
        for axis,key,title in [(axes[0,0],'rmsRho','Density residual change'),(axes[0,1],'rmsRhoE','Energy residual change')]:
            axis.plot(h['Inner_Iter'],h[key]-h[key][0],label=name,color=color,lw=1);axis.set_title(title)
            axis.set_xlabel('Local iteration');axis.set_ylabel('log10(R / R at iteration 0)')
        xx=np.unique(v[:,0]);xx=xx[(xx>0)&(xx<max(xx))];s=(.446278-xx)/.0254
        for axis,col,title,unit in [(axes[1,0],0,'Circumferential mean wall pressure','p (Pa)'),(axes[1,1],2,'Wall-fluid temperature, adiabatic','T (K)')]:
            mean=[np.mean(nodal[v[:,0]==x,col]) for x in xx]
            axis.plot(s,mean,label=name,color=color,lw=1);axis.set_title(title);axis.set_xlabel('Axial distance from nose (in)');axis.set_ylabel(unit)
    for axis in axes.flat:axis.grid(alpha=.2);axis.legend(frameon=False,fontsize=8)
    fig.suptitle('Unconverged numerical study - not experimental validation',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(str(ROOT/'docs/validation/v116/collar_comparison.png'),dpi=160);plt.close(fig)
    print(json.dumps(dict(cases=records,difference=difference,nominalStationMatchedRecords=matches)))
if __name__=='__main__':main()
