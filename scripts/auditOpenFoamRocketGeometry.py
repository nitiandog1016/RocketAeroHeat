"""独立核对原全壁及原生失败指标，定位具体单元；只读、不生成优化候选。"""
from pathlib import Path
import json,re,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from openFoamConservation import geometry
from openFoamFields import internal_field
from openFoamAirThermo import evaluate


def internal_area_determinant(sf,owner,neighbor,count):
    """模块一：原生三维内部/耦合面面积张量的行列式，普通边界面不参与。"""
    n=len(neighbor);a=sf[:n];mag=np.linalg.norm(a,axis=1)
    counts=np.bincount(np.r_[owner[:n],neighbor],minlength=count)
    areas=np.bincount(np.r_[owner[:n],neighbor],weights=np.r_[mag,mag],minlength=count)
    mean=np.divide(areas,counts,out=np.zeros(count),where=counts>0)
    tensor=np.zeros((count,3,3));outer=a[:,:,None]*a[:,None,:]
    np.add.at(tensor,owner[:n],outer);np.add.at(tensor,neighbor,outer)
    tensor=np.divide(tensor,mean[:,None,None]**2,out=np.zeros_like(tensor),where=mean[:,None,None]>0)
    return abs(np.linalg.det(tensor)),counts


def main():
    start=time.monotonic();base=ROOT/'outputs/v200_validation';out=base/'rocket_geometry01';out.mkdir(exist_ok=False);case=base/'rocket02'
    shutil.copyfile(__file__,out/'executed_entry.py')
    source=ROOT/'outputs/v198_validation/prepared_audit01/primitiveMeshCheck.C';shutil.copyfile(source,out/source.name)
    assert sha(source)==sha(ROOT/'docs/validation/v198/evidence/prepared_audit01/primitiveMeshCheck.C')
    g=geometry(case);p=g['points'];c=g['cells'];ni=len(g['neighbor']);metric,counts=internal_area_determinant(g['sf'],g['owner'],g['neighbor'],len(c))
    check=json.loads((base/'rocket_check02/report.json').read_text());assert not check['passed'] and check['process']['passed']
    log=(base/'rocket_check02/native/stdout.txt').read_text();native_min=float(re.search(r'wellposedness\) : minimum: ([\deE.+-]+)',log).group(1))
    error=abs(metric.min()/native_min-1);bad=np.flatnonzero(metric<.001)
    assert len(bad)==19 and error<2e-10
    surface=json.loads((ROOT/'outputs/v194_validation/rocket_design03/surface.json').read_text());original=np.array(surface['facesOneBased'])-1
    wall=next(q for q in g['patches'] if q['name']=='wall');ids=np.arange(wall['startFace'],wall['startFace']+wall['nFaces'])
    assert np.array_equal(p[:9878],surface['vertices']) and np.array_equal(np.sort(original,axis=1),np.sort(g['faces'][ids],axis=1))
    xyz=p[original];edges=np.stack([np.linalg.norm(xyz[:,1]-xyz[:,0],axis=1),np.linalg.norm(xyz[:,2]-xyz[:,1],axis=1),np.linalg.norm(xyz[:,0]-xyz[:,2],axis=1)],axis=1)
    semi=edges.sum(axis=1)/2;areas=np.sqrt(semi*np.prod(semi[:,None]-edges,axis=1))
    area_error=float(abs(areas.sum()-1.9646292565715553));assert area_error<2e-13
    v=p[c];vol=abs(np.einsum('ij,ij->i',np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),v[:,3]-v[:,0]))/6
    volume_error=float(abs(vol.sum()/381.4346146072696-1));assert volume_error<2e-13
    # 位置由原体积分点给出，原壁点/外界点固定与内部自由点分别登记。
    boundary_nodes=set(g['faces'][ni:].reshape(-1).tolist());wall_nodes=set(original.reshape(-1).tolist());rows=[]
    for cell in bad:
        nodes=c[cell];touch=ids[g['owner'][ids]==cell]
        rows.append(dict(cellZeroBased=int(cell),metric=float(metric[cell]),centerM=g['centers'][cell].tolist(),
            nodesZeroBased=nodes.tolist(),internalFaces=int(counts[cell]),wallVertexCount=sum(int(n) in wall_nodes for n in nodes),
            fixedBoundaryVertexCount=sum(int(n) in boundary_nodes for n in nodes),
            originalWallFaceIDs=(touch-wall['startFace']+1).tolist()))
    f={name:internal_field(case/'0'/name,len(c),name=='U') for name in ['p','T','U']};gas=evaluate(f['p'],f['T'])
    mach=np.linalg.norm(f['U'],axis=1)/gas['soundSpeedMps'];assert np.max(abs(mach-5))<2e-13
    report=dict(passed=True,scope='independent geometry and failed native metric; not mesh acceptance',
        nativeMeshAccepted=False,originalWallExact=True,originalWallFaces=len(original),originalAreaM2=float(areas.sum()),
        areaAbsoluteErrorM2=area_error,volumeRelativeError=volume_error,nativeDeterminantMin=native_min,
        replayDeterminantMin=float(metric.min()),determinantMinRelativeError=float(error),badCells=rows,
        initialMachMaxError=float(np.max(abs(mach-5))),minimumInternalFaces=int(counts.min()),
        newQualifiedHeatAreaM2=0,originalRocketCFDRuns=0,sourceSHA256={str(q):sha(q) for q in [source,case/'exact_mesh.npz',base/'rocket_check02/report.json']},
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','areaAbsoluteErrorM2','volumeRelativeError','determinantMinRelativeError','minimumInternalFaces','elapsedSeconds']}))
    print(json.dumps(dict(badCells=len(rows),wallVertices=[r['wallVertexCount'] for r in rows],fixedBoundaryVertices=[r['fixedBoundaryVertexCount'] for r in rows],xRange=[min(r['centerM'][0] for r in rows),max(r['centerM'][0] for r in rows)])))


if __name__=='__main__':main()
