"""只读定位原层与原生门槛的不相容，并仅修复格式排序；不改物理网格。"""
import json,time,shutil,re
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha,header
from auditOpenFoamRocketGeometry import internal_area_determinant
BASE=ROOT/'outputs/v208_validation'

def main():
    start=time.monotonic();out=BASE/'layer_audit01';out.mkdir(exist_ok=False)
    case=BASE/'layer_mesh01';d=np.load(case/'exact_mesh.npz');p=d['points'];c=d['cells'];faces=d['faces'];size=d['faceSizes'];owner=d['owner'];nei=d['neighbour'];ni=len(nei)
    sf=np.zeros((len(faces),3))
    for n in [3,4]:
        mask=size==n;q=p[faces[mask,:n]];sf[mask]=np.cross(q,np.roll(q,-1,axis=1)).sum(axis=1)/2
    det,internal_count=internal_area_determinant(sf,owner,nei,len(c));bad=det<.001
    log=(BASE/'layer_check01/native/stdout.txt').read_text()
    native_min=float(re.search(r'wellposedness\) : minimum: ([\deE.+-]+)',log).group(1))
    native_bad=int(re.search(r'Cells with small determinant.*number of cells: (\d+)',log).group(1))
    assert int(bad.sum())==native_bad==139027 and abs(det.min()/native_min-1)<5e-6
    # 模块一：内部面只进行编号置换；原壁面的索引和对应结果位置不变。
    order=np.r_[np.lexsort((nei,owner[:ni])),np.arange(ni,len(owner))]
    new_faces=faces[order];new_owner=owner[order];new_nei=nei[order[:ni]];new_sizes=size[order]
    assert np.array_equal(order[ni:],np.arange(ni,len(owner)))
    assert np.all((new_owner[1:ni]>new_owner[:ni-1])|((new_owner[1:ni]==new_owner[:ni-1])&(new_nei[1:]>new_nei[:-1])))
    det2,_=internal_area_determinant(sf[order],new_owner,new_nei,len(c))
    assert np.array_equal(bad,det2<.001) and np.max(abs(det2-det)/np.maximum(det,1e-30))<2e-9
    destination=BASE/'layer_mesh02';mesh=destination/'constant/polyMesh';mesh.mkdir(parents=True);(destination/'system').mkdir();(destination/'0').mkdir()
    # 几何只读复制到新命名目录，原记录完整保留，避免重跑转换/原生诊断。
    for name in ['points','boundary']:shutil.copyfile(case/'constant/polyMesh'/name,mesh/name)
    for name in ['controlDict','fvSchemes','fvSolution']:shutil.copyfile(case/'system'/name,destination/'system'/name)
    def listing(name,kind,rows,count):
        with (mesh/name).open('w',encoding='ascii',newline='\n') as f:
            f.write(header(name,kind)+str(count)+'\n(\n')
            for row in rows:f.write(row+'\n')
            f.write(')\n// ************************************************************************* //\n')
    listing('faces','faceList',(str(n)+'('+' '.join(map(str,f[:n]))+')' for f,n in zip(new_faces,new_sizes)),len(new_faces))
    listing('owner','labelList',map(str,new_owner),len(new_owner));listing('neighbour','labelList',map(str,new_nei),len(new_nei))
    kwargs={k:d[k] for k in d.files};kwargs.update(faces=new_faces,faceSizes=new_sizes,owner=new_owner,neighbour=new_nei)
    np.savez_compressed(destination/'exact_mesh.npz',**kwargs)
    # 模块二：原面ID/节点/Heron面积独立对照，按实际层定位失败，不重释为分离区。
    source=json.loads((ROOT/'outputs/v194_validation/rocket_design03/surface.json').read_text());v=np.array(source['vertices']);f=np.array(source['facesOneBased'])-1
    ids=d['originalWallFaceIDs'];original=d['originalNodeIDs'];wall=faces[ni:ni+7838,:3]
    assert np.array_equal(np.sort(original[wall]-1,axis=1),np.sort(f[ids-1],axis=1))
    assert np.array_equal(p[wall],v[original[wall]-1])
    xyz=p[wall];edges=np.linalg.norm(xyz-np.roll(xyz,1,axis=1),axis=2);sp=edges.sum(axis=1)/2
    area=np.sqrt(sp*np.prod(sp[:,None]-edges,axis=1));assert abs(area.sum()-.7198603834559664)<1e-13
    layer=((original[c]-1)//9878).min(axis=1);assert np.array_equal(np.bincount(layer),np.full(32,7838))
    pure=internal_count==5
    records=[dict(layer=int(i),cells=int((layer==i).sum()),determinantFailed=int((bad&(layer==i)).sum()),
                  fullyInternalFailed=int((bad&pure&(layer==i)).sum())) for i in range(32)]
    wall_info=[]
    for j,faceid in enumerate(ids):
        cell=owner[ni+j];wall_info.append([int(faceid),float(det[cell]),int(bad[cell])])
    np.savetxt(out/'wall_adjacent_determinant.csv',np.array(wall_info),delimiter=',',header='originalFaceID,determinant,failed',comments='',fmt=['%d','%.17g','%d'])
    source_path=ROOT/'outputs/v198_validation/prepared_audit01/primitiveMeshCheck.C'
    shutil.copyfile(source_path,out/'primitiveMeshCheck.C');shutil.copyfile(__file__,out/'executed_entry.py')
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert all(sha(ROOT/n.replace('\\','/'))==h for n,h in protected.items())
    result=dict(passed=True,scope='ordering correction and native determinant replay; not native mesh acceptance',
        originalWallFaces=7838,originalWallAreaM2=float(area.sum()),layers=records,minimumDeterminant=float(det.min()),
        nativePrintedMinimum=native_min,printedMinimumRelativeDifference=float(abs(det.min()/native_min-1)),
        failedDeterminantCells=int(bad.sum()),fullyInternalFailedCells=int((bad&pure).sum()),
        originalNativeMeshAccepted=False,correctedMeshNativeCheckPerformed=False,geometryQualityAccepted=False,
        orderingCorrected=True,allBoundaryFaceIndicesPreserved=True,originalCellsAndPointsPreserved=True,
        determinantPermutationMaxRelative=float(np.max(abs(det2-det)/np.maximum(det,1e-30))),
        sourceSHA256={str(source_path.relative_to(ROOT)):sha(source_path),'layer_mesh01/exact_mesh.npz':sha(case/'exact_mesh.npz')},
        correctedSHA256={p.name:sha(p) for p in mesh.iterdir()},protectedFilesUnchanged=29,
        nativeRuns=0,CFDRuns=0,physicalQualification=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    assert result['elapsedSeconds']<60
    save(destination/'report.json',result);save(out/'report.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['layers','sourceSHA256','correctedSHA256']}))

if __name__=='__main__':main()
