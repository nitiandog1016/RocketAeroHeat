"""模块职责：独立读回真实区域，核对拓扑/双侧特征/源状态，保留元数据错误。"""
from pathlib import Path
import sys,struct
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq
from runSU2RocketPhaseV222 import ROOT,read,save,sha
import prepareViscousRegionV244 as run
BASE=run.BASE;DOC=run.DOC;CLOSED=True


def parse_mesh(path):
    """独立文本读取器：保留每种单元/面型、局部号与顺序，不借用生产数组。"""
    with Path(path).open(encoding='utf8') as f:
        assert next(f).strip()=='NDIME= 3'
        n=int(next(f).split('=')[1]);cells=np.full((n,6),-1,dtype=np.int64);types=[]
        for i in range(n):
            row=list(map(int,next(f).split()));assert row[-1]==i and row[0] in (10,13)
            assert len(row)==(6 if row[0]==10 else 8)
            types.append(row[0]);cells[i,:len(row)-2]=row[1:-1]
        n=int(next(f).split('=')[1]);xyz=np.empty((n,3))
        for i in range(n):
            row=next(f).split();assert int(row[3])==i;xyz[i]=list(map(float,row[:3]))
        nm=int(next(f).split('=')[1]);markers={}
        for j in range(nm):
            name=next(f).split('=')[1].strip();n=int(next(f).split('=')[1]);faces=np.full((n,4),-1,dtype=np.int64)
            for k in range(n):
                row=list(map(int,next(f).split()));assert row[0] in (5,9)
                assert len(row)==(4 if row[0]==5 else 5);faces[k,:len(row)-1]=row[1:]
            assert name not in markers;markers[name]=faces
        assert f.read().strip()==''
    return xyz,cells,np.array(types),markers


def audit():
    """模块一：原号、所有区域体单元和完整边界从文本独立恢复。"""
    out=BASE/'audit02';out.mkdir(exist_ok=False);src=BASE/'region02'
    result=read(src/'report.json');assert result['passed']
    xyz,cells,types,markers=parse_mesh(src/'region_geometry.su2');z=np.load(str(src/'region.npz'))
    assert np.array_equal(xyz,z['vertices']) and np.array_equal(cells,z['localCells'])
    d=loadmat(str(run.VOLUME),variable_names=['vertices','wedges','tetrahedra'])
    ids=z['globalNodeIDs'];cid=z['motherCellIDs'];nw=len(d['wedges'])
    assert np.array_equal(xyz,d['vertices'][ids]) and np.all(types==np.where(cid<nw,13,10))
    true=np.full_like(cells,-1);prism=cid<nw
    true[prism]=d['wedges'][cid[prism]]-1;true[~prism,:4]=d['tetrahedra'][cid[~prism]-nw]-1
    recovered=cells.copy();valid=recovered>=0;recovered[valid]=ids[recovered[valid]]
    assert np.array_equal(recovered,true)
    half=[]
    for typ,patterns in ((13,((0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5))),
                         (10,((0,2,1),(0,1,3),(1,2,3),(2,0,3)))):
        q=cells[types==typ]
        for pat in patterns:
            a=np.full((len(q),4),-1,dtype=np.int64);a[:,:len(pat)]=q[:,pat];half.append(a)
    keys,count=np.unique(np.sort(np.concatenate(half),axis=1),axis=0,return_counts=True)
    assert count.max()==2
    boundary=np.sort(np.concatenate(list(markers.values())),axis=1);boundary=np.unique(boundary,axis=0)
    assert len(boundary)==sum(len(v) for v in markers.values())
    assert np.array_equal(boundary,keys[count==1])
    for name,face in markers.items():assert np.array_equal(face,z[name+'LocalFaces'])
    # 独立读取二进制，按全局号取得原状态；研究场不当物理解。
    with run.FRAME.open('rb') as f:
        magic,n,k=struct.unpack('<8sQQ',f.read(24));assert magic==b'RAH235LE' and (n,k)==(335726,34)
        original=np.fromfile(f,dtype='<f8').reshape(n,k)
    original=original[np.argsort(original[:,0])];state=z['regionState']
    assert np.array_equal(state.view(np.uint64),original[ids].copy().view(np.uint64))
    """模块二：由表数据直接积分Cv，Brent反解，不调用生产的二次根反解器。"""
    gas=read(ROOT/'data/thermo/phase_aware_dry_air_v1.json');T=np.array(gas['temperatureK'])
    cp=np.array(gas['cpJkgK']);R=gas['gasConstantJkgK'];cv=cp-R
    knots=np.r_[0,np.cumsum(np.diff(T)*(cv[:-1]+cv[1:])/2)]
    ref=int(np.flatnonzero(T==298.15)[0]);knots-=knots[ref]+R*298.15
    def energy(temp):
        i=min(len(T)-2,int(np.searchsorted(T,temp,side='right')-1));dt=temp-T[i]
        return knots[i]+cv[i]*dt+.5*(cv[i+1]-cv[i])/(T[i+1]-T[i])*dt*dt
    table=np.loadtxt(str(src/'interface.csv'),delimiter=',',skiprows=1)
    idx={k:v for k,v in np.load(str(run.INDEX)).items()}
    # npz按键访问会重新解压；大母面表只读一次，禁止在每面循环重复解压。
    mother_faces=np.load(str(run.PART))['motherFaces']
    maxerror=0.;counts=[];wall=set(markers['wall'].ravel())-{-1};contact=[]
    for j,face in enumerate(markers['unclosed_internal_interface']):
        nodes=face[face>=0];point=xyz[nodes];centroid=point.mean(axis=0)
        # 将多边形扇形三角化；与生产的多边形环叉积独立对照。
        normal=sum(np.cross(point[i]-point[0],point[i+1]-point[0])/2 for i in range(1,len(point)-1))
        up=int(table[j,1]);cc=(d['wedges'][up]-1) if up<nw else (d['tetrahedra'][up-nw]-1)
        if np.dot(normal,centroid-d['vertices'][cc].mean(axis=0))<0:normal=-normal
        area=np.linalg.norm(normal);normal/=area
        u=state[nodes,5:10].mean(axis=0);v=u[1:4]/u[0];target=u[4]/u[0]-.5*v@v-300000
        temp=brentq(lambda tt:energy(tt)-target,T[0],T[-1],xtol=1e-10)
        pressure=u[0]*R*temp;sound=np.sqrt(np.interp(temp,T,cp)/np.interp(temp,T,cv)*R*temp)
        mn=float(v@normal/sound);incoming=0 if mn>=1 else (1 if mn>0 else (4 if mn>-1 else 5))
        counts.append(incoming);contact.append(bool(set(nodes)&wall))
        expected=np.r_[centroid,normal,area,pressure,temp,mn]
        err=float(max(abs(table[j,4:14]-expected)/np.maximum(1,abs(expected))))
        maxerror=max(maxerror,err);assert err<2e-9
        assert incoming==int(table[j,14]) and table[j,14]+table[j,15]==5 and table[j,16]==0
        assert int(table[j,0])==idx['interfaceMotherFaces'][j]
        assert int(table[j,1])==idx['upstreamInterfaceCells'][j] and int(table[j,2])==idx['aftInterfaceCells'][j]
        mother=mother_faces[int(table[j,0])]
        assert np.array_equal(np.sort(ids[nodes]),np.sort(mother[mother>=0]))
    contact=int(sum(contact));assert contact==46
    # 元数据修复只读旧场；首报告14等于严格零壁距数，不能追改成拓扑接触数。
    corrected=run.wall_contact_stats(markers['unclosed_internal_interface'],markers['wall'],table)
    assert corrected['interfaceFacesTouchingWall']==contact
    fixed=dict(result);fixed['primitiveFaceCharacteristics']=dict(result['primitiveFaceCharacteristics'],**corrected)
    fixed['metadataCorrection']=dict(originalReportSHA256=sha(src/'report.json'),
        originalWallContactCount=result['primitiveFaceCharacteristics']['interfaceFacesTouchingWall'],
        correctedWallContactCount=contact,originalReportPreserved=True,noNewFlow=True)
    save(out/'corrected_region_report.json',fixed)
    report=dict(passed=True,allRegionalCellsTextRoundtrip=True,allRegionalCoordinatesBitEqual=True,
        originalStateScalarsBitEqual=int(state.size),boundaryCompleteNoDuplicates=True,
        interfaceFaces=1984,independentThermodynamicAndNormalMaxScaledError=maxerror,
        incomingUpstreamCounts={str(i):counts.count(i) for i in (0,1,4,5)},
        exactTopologicalWallContactFaces=contact,strictlyZeroStoredDistanceFaces=corrected['interfaceFacesWithExactlyZeroStoredWallDistance'],
        maximumStoredDistanceAtContactM=float(max(table[np.isin(markers['unclosed_internal_interface'],list(wall)).any(axis=1),17])),
        actualCorrectedMetadataFunctionTested=True,sourceFieldChanged=False,
        productionAdmitted=False,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'report.json',report);print(report)


def main():
    if CLOSED:raise SystemExit('v244 region audit closed; use NEXT_DECISION')
    run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC;size=run.guard.size_bytes
    run.guard.size_bytes=lambda p:size(BASE) if Path(p)==BASE or Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
        run.VOLUME,run.FRAME,run.INDEX,run.PART,ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    locks += [p for p in (BASE/'region02').iterdir() if p.is_file()]
    return run.guard.supervise('audit2',[sys.executable,str(Path(__file__).resolve()),'task'],60,256,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v244 region audit closed; use NEXT_DECISION')
    if sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
