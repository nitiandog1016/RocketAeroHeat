"""模块职责：原物理九点入口不变，比较新网格真实场的流管、热状态和横向度量。"""
from pathlib import Path
import sys
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
import runSU2RocketPhaseV222 as guard
from tetraStreamtubes import TetraField,tube,LABELS
from phaseEulerRegionalVolume import read_tetra_mesh
from auditForebodyStreamtubesV241 import independent_states,segment_boundary_hits
BASE=ROOT/'outputs/v242_validation';DOC=ROOT/'docs/validation/v242';CASE=BASE/'run01'
OLD=ROOT/'outputs/v241_validation/streamtubes02';SOURCE=BASE/'tubes01'
CLOSED=True


def produce():
    """模块一：逐位复用已记录的实际入口，不按新网格质心移动种子。"""
    audit=read(BASE/'flow_audit01/report.json');assert audit['passed'] and audit['completedUpdates']==400
    original=read(OLD/'manifest.json');x,t,m=read_tetra_mesh(CASE/'mesh.su2')
    a=np.loadtxt(str(CASE/'state_000400.csv'),delimiter=',',skiprows=1);a=a[np.argsort(a[:,0])]
    assert np.array_equal(a[:,1:4],x)
    field=TetraField(x,t,a[:,4:9]);SOURCE.mkdir(exist_ok=False)
    manifest=dict(passed=False,status='unqualified-research',physicalQualification=False,edgeAccepted=False,heatProviderCalled=False,
      newQualifiedAreaM2=0,sourceUpdate=400,sourceBranch='phase-air-single-shell-forebody-euler-v242',
      sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [CASE/'mesh.su2',CASE/'state_000400.csv',CASE/'original_faces.csv',
        OLD/'manifest.json',BASE/'flow_audit01/report.json',ROOT/'data/thermo/phase_aware_dry_air_v1.json']},
      dsM=original['dsM'],requestedSteps=original['requestedSteps'],samePhysicalNineSeeds=True,
      originalFaces=19752,originalAreaM2=1.9646292565715553,newStepConvergenceChecked=False,tubes=[])
    for info in original['tubes']:
        seed=np.array(info['seedXYZ']);beta=np.array(info['seedLabelBeta']);eta=np.array(info['seedLabelEta']);delta=info['labelDeltaM']
        owner,_=field.locate(seed,None);rays,metric=tube(field,seed,beta,eta,delta,manifest['dsM'],manifest['requestedSteps'])
        arrays={}
        for j,r in enumerate(rays):
            assert np.array_equal(r['xyz'][0],seed+delta*(LABELS[j,0]*beta+LABELS[j,1]*eta))
            for name in ('xyz','cell','weights','U','s'):arrays['ray%d_%s'%(j,name)]=r[name]
        for name,value in metric.items():arrays['metric_'+name]=np.asarray(value)
        path=SOURCE/info['file'];np.savez_compressed(str(path),**arrays)
        record={key:info[key] for key in ('id','originalSeedFaceID','seedXYZ','seedLabelBeta','seedLabelEta','labelDeltaM','file')}
        record.update(newOwnerCell=owner,commonStations=metric['commonStations'],rayStationCounts=[len(r['s']) for r in rays],
          rayStops=[{key:r[key] for key in ('stopReason','stopPhase','stopPoint','stopBoundaryNodes')} for r in rays],
          maxAbsHRelativeDrift=float(max(abs(metric['HRelativeDrift']))),maxAbsEntropyDriftJkgK=float(max(abs(metric['entropyDrift']))),
          massFluxRatioRange=[float(min(metric['massFluxRatio'])),float(max(metric['massFluxRatio']))],
          validMetricStations=int(sum(metric['valid'])),sha256=sha(path))
        manifest['tubes'].append(record);save(SOURCE/'manifest.json',manifest);print('tube',info['id'],'common',metric['commonStations'],flush=True)
    manifest['passed']=True;save(SOURCE/'manifest.json',manifest)


def audit_export():
    """模块二：独立审核全部保存站、接受线段与边界交点，再导出对照CSV/科研图。"""
    original=read(OLD/'manifest.json');manifest=read(SOURCE/'manifest.json');assert manifest['passed']
    for p,h in manifest['sourceSHA256'].items():assert sha(ROOT/p)==h
    x,t,markers=read_tetra_mesh(CASE/'mesh.su2');raw=np.loadtxt(str(CASE/'state_000400.csv'),delimiter=',',skiprows=1);raw=raw[np.argsort(raw[:,0])]
    lookup={tuple(sorted(f)):tag for tag,faces in markers.items() for f in faces};records=[];count=0;maxu=maxm=0.;stops=[]
    out=BASE/'tube_audit01';out.mkdir(exist_ok=False);header='sM,xM,yM,zM,pPa,TK,HJkg,SJkgK,hBeta,hEta,massFluxRatio,HRelativeDrift,entropyDrift,metricValid'
    for info,old in zip(manifest['tubes'],original['tubes']):
        for key in ('seedXYZ','seedLabelBeta','seedLabelEta','labelDeltaM','originalSeedFaceID'):assert info[key]==old[key]
        assert sha(SOURCE/info['file'])==info['sha256'];z=np.load(str(SOURCE/info['file']));prior=np.load(str(OLD/old['file']))
        n=info['commonStations'];states=[]
        for j in range(9):
            p=z['ray%d_xyz'%j];cells=z['ray%d_cell'%j];u=z['ray%d_U'%j];q=x[t[cells]]
            assert np.array_equal(p[0],prior['ray%d_xyz'%j][0]) and len(p)==info['rayStationCounts'][j]
            v=np.linalg.solve((q[:,1:]-q[:,0,None]).transpose(0,2,1),p-q[:,0]);w=np.column_stack([1-v.sum(axis=1),v])
            assert np.max(abs(w-z['ray%d_weights'%j]))<5e-10 and w.min()>=-2e-11
            expected=np.einsum('ni,nij->nj',w,raw[t[cells],4:9]);err=float(np.max(abs(expected-u)/np.maximum(1,abs(expected))))
            assert err<2e-10;maxu=max(maxu,err);count+=len(p);states.append(independent_states(u))
            assert np.array_equal(z['ray%d_s'%j],np.arange(len(p))*manifest['dsM'])
            assert not segment_boundary_hits(p,x,markers),(info['id'],j)
            stop=info['rayStops'][j]
            if stop['stopReason']=='TUBE_BOUNDARY_CROSSING':
                face=stop['stopBoundaryNodes'];tag=lookup[tuple(sorted(face))];pt=np.asarray(stop['stopPoint']);tri=x[face]
                b=np.linalg.lstsq((tri[1:]-tri[0]).T,pt-tri[0],rcond=None)[0]
                assert min(b)>=-1e-10 and sum(b)<=1+1e-10 and np.linalg.norm(pt-(tri[0]+b@(tri[1:]-tri[0])))<1e-11
                stops.append(dict(tube=info['id'],ray=j,marker=tag,xyz=pt.tolist()))
            else:assert stop['stopReason']=='length-limit'
        center=states[0];delta=info['labelDeltaM'];direction=center['v'][:n]/center['U'][:n,None]
        jb=(z['ray3_xyz'][:n]-z['ray2_xyz'][:n])/delta;je=(z['ray7_xyz'][:n]-z['ray6_xyz'][:n])/delta
        area=np.array([np.linalg.det(np.vstack([a,b,c])) for a,b,c in zip(direction,jb,je)])
        flux=center['rho'][:n]*center['U'][:n]*area
        for key,value in [('hBeta',np.linalg.norm(jb-direction*np.sum(jb*direction,axis=1)[:,None],axis=1)),
          ('signedArea',area),('massFluxRatio',flux/flux[0]),('HRelativeDrift',(center['H'][:n]-center['H'][0])/abs(center['H'][0])),
          ('entropyDrift',center['S'][:n]-center['S'][0])]:
            err=float(np.max(abs(value-z['metric_'+key])/np.maximum(1,abs(value))));assert err<2e-10;maxm=max(maxm,err)
        data=np.column_stack([z['metric_s'],z['ray0_xyz'][:n],center['p'][:n],center['T'][:n],center['H'][:n],center['S'][:n]]+
          [z['metric_'+key] for key in ('hBeta','hEta','massFluxRatio','HRelativeDrift','entropyDrift','valid')])
        path=out/('tube%02d.csv'%info['id']);np.savetxt(str(path),data,delimiter=',',header=header,comments='',fmt='%.17g')
        assert np.array_equal(np.loadtxt(str(path),delimiter=',',skiprows=1),data)
        overlap=min(n,old['commonStations']);records.append(dict(id=info['id'],commonStations=n,oldCommonStations=old['commonStations'],comparisonPrefixStations=overlap,
          maxAbsHRelativeDrift=info['maxAbsHRelativeDrift'],oldMaxAbsHRelativeDrift=old['maxAbsHRelativeDrift'],
          massFluxRatioRange=info['massFluxRatioRange'],oldMassFluxRatioRange=old['massFluxRatioRange'],
          samePrefixMaxHDriftNew=float(max(abs(z['metric_HRelativeDrift'][:overlap]))),
          samePrefixMaxHDriftOld=float(max(abs(prior['metric_HRelativeDrift'][:overlap]))),csvValues=int(data.size)))
    # 模块三：相同种子位置的新旧沿程对照，不用不同管长的极值冒充同区间精度比较。
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7));fig.subplots_adjust(left=.09,right=.97,bottom=.12,top=.86,wspace=.28,hspace=.35)
    for ax,info,record in zip(axs.ravel(),manifest['tubes'],records):
        for label,folder,style,color in [('old mesh / update 200',OLD,'--','#7c858c'),('new mesh / update 400',SOURCE,'-','#206987')]:
            z=np.load(str(folder/info['file']));ax.plot(z['metric_s'],100*z['metric_HRelativeDrift'],style,color=color,lw=1.3,label=label)
        ax.set_title('Tube %d | original face %d'%(info['id'],info['originalSeedFaceID']),loc='left',fontsize=10)
        ax.set_xlabel('Arc length [m]');ax.set_ylabel('Total H drift from seed [%]');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Same physical nine-ray seeds | new saved Euler field comparison\nResearch only: near-wall P1 streamtubes are NOT accepted boundary-layer edges',fontsize=12)
    fig.text(.09,.025,'Solid / dashed endpoints are actual common-prefix stops. No seed shift, projection or filled continuation.\nChanges include both the grid and achieved iterative convergence; this is not a grid-independence study.',fontsize=9,color='#444')
    fig.savefig(str(out/'tube_comparison.png'),dpi=160);plt.close(fig)
    save(out/'report.json',dict(passed=True,savedRayStationsAudited=count,maxIndependentFiveStateScaledDifference=maxu,
      maxIndependentMetricScaledDifference=maxm,records=records,boundaryStops=stops,samePhysicalNineSeeds=True,
      allSavedNineRaySegmentsChecked=True,exactCSVRoundtripValues=sum(r['csvValues'] for r in records),
      totalCommonStations=sum(r['commonStations'] for r in records),newStepConvergenceChecked=False,
      physicalQualification=False,edgeAccepted=False,newQualifiedAreaM2=0))


def main(stage):
    if CLOSED:raise SystemExit('v242 tube comparison closed; use NEXT_DECISION')
    guard.BASE=BASE;guard.CASE=CASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',DOC/'STREAMTUBE_PLAN.md',ROOT/'scripts/tetraStreamtubes.py',
      ROOT/'scripts/auditForebodyStreamtubesV241.py',ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/runSU2RocketPhaseV222.py',
      ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json',
      BASE/'flow_audit01/report.json']+[CASE/p for p in ('mesh.su2','state_000400.csv','original_faces.csv')]+list(OLD.glob('*'))
    if stage=='audit':locks+=list(SOURCE.glob('*'))
    assert stage in ('produce','audit')
    return guard.supervise('tubes_'+stage,[sys.executable,str(Path(__file__).resolve()),'task_'+stage],60,256,locks,BASE,'offline' if stage=='produce' else 'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v242 tube comparison closed; use NEXT_DECISION')
    if sys.argv[1]=='task_produce':produce()
    elif sys.argv[1]=='task_audit':audit_export()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
