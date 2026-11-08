"""只读验证v207实际热壁输出和原壁区域；损坏样本只存在内存中。"""
from pathlib import Path
import io,json,re,shutil,sys,time
from unittest.mock import patch
from collections import Counter
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepareOpenFoamControl import sha,save
from readOpenFoamWallHeat import read_wall_heat

BASE=ROOT/'outputs/v207_validation'
CASE=BASE/'thermal03'

def reader_checks():
    """模块一：对真实保存帧核对符号/身份，注入损坏后必须明确拒绝。"""
    old=json.loads((BASE/'thermal_audit03/report.json').read_text())
    for frame in old['records']:
        result=read_wall_heat(CASE,frame['time'],32)
        assert result['rows']==frame['fields']['rows']
        rows=result['rows']
        assert rows[0]['heatFluxIntoWallWm2']>0 and rows[1]['heatFluxIntoWallWm2']<0
        assert [r['ownerCellIndex'] for r in rows]==[0,31]
        assert [r['nativeFaceIndex'] for r in rows]==[31,32]
        assert not result['physicalQualification'] and not result['originalRocketMapped']
    actual_read=Path.read_text;failures=[]
    def rejected(label,relative,transform):
        target=CASE/relative
        def injected(path,*args,**kwargs):
            text=actual_read(path,*args,**kwargs)
            return transform(text) if path==target else text
        with patch.object(Path,'read_text',injected):
            try:read_wall_heat(CASE,'4e-08',32)
            except (ValueError,FileNotFoundError) as error:failures.append(dict(name=label,error=str(error)))
            else:raise AssertionError('accepted damaged input: '+label)
    def replace_once(pattern,value):
        def change(text):
            changed,n=re.subn(pattern,value,text,count=1)
            assert n==1 and changed!=text
            return changed
        return change
    for field in ['p','T','wallHeatFlux']:
        rel='4e-08/'+field
        rejected(field+' truncated',rel,lambda t:t[:-90])
        rejected(field+' dimensions',rel,replace_once(r'dimensions\s+\[[^\]]+\]','dimensions [0 0 0 0 0 0 0]'))
        rejected(field+' object',rel,replace_once(r'object\s+\w+;','object wrong;'))
        rejected(field+' time',rel,replace_once(r'location\s+"4e-08";','location "3e-08";'))
        rejected(field+' class',rel,replace_once(r'class\s+volScalarField;','class volVectorField;'))
    rejected('nonfinite heat','4e-08/wallHeatFlux',replace_once(r'value\s+uniform\s+[^;]+;','value uniform nan;'))
    rejected('heat count','4e-08/wallHeatFlux',replace_once(r'value\s+uniform\s+[^;]+;','value nonuniform List<scalar> 2 (1 2);'))
    rejected('missing heat value','4e-08/wallHeatFlux',replace_once(r'value\s+uniform\s+[^;]+;',''))
    rejected('wall temperature identity','4e-08/T',replace_once(r'type\s+fixedValue;','type zeroGradient;'))
    rejected('wall temperature domain','4e-08/T',replace_once(r'value\s+uniform\s+300;','value uniform 50;'))
    rejected('wall temperature negative','4e-08/T',replace_once(r'value\s+uniform\s+300;','value uniform -1;'))
    rejected('wall overlap','constant/polyMesh/boundary',replace_once(r'startFace\s+32;','startFace 31;'))
    rejected('wall past end','constant/polyMesh/boundary',replace_once(r'startFace\s+32;','startFace 99999;'))
    rejected('owner count','constant/polyMesh/owner',replace_once(r'\n161\s*\(','\n162\n('))
    rejected('owner range','constant/polyMesh/owner',replace_once(r'\n161\s*\(\s*0','\n161\n(\n32'))
    rejected('missing wall','constant/polyMesh/boundary',lambda s:re.sub(r'type\s+wall;','type patch;',s))
    rejected('wrong caloric table','constant/physicalProperties',replace_once('1002.2432415219383','1003.2432415219383'))
    rejected('wrong energy identity','constant/physicalProperties',replace_once('energy rahPhaseAwareInternalEnergy','energy sensibleInternalEnergy'))
    rejected('unvalidated turbulence','constant/momentumTransport',replace_once('simulationType laminar','simulationType RAS'))
    rejected('unvalidated conduction','constant/thermophysicalTransport',replace_once('model Fourier','model unknown'))
    rejected('wrong sign contract','thermo_contract.json',lambda s:s.replace('"nativeWallHeatFluxToIntoWallSign": -1','"nativeWallHeatFluxToIntoWallSign": 1'))
    for bad in ['0','nan','../4e-08','-1']:
        try:read_wall_heat(CASE,bad,32)
        except ValueError as error:failures.append(dict(name='bad time '+bad,error=str(error)))
        else:raise AssertionError('accepted bad time')
    # 有限零热流是有效物理数值，不能与缺失值混淆；这里只改内存样本。
    target=CASE/'4e-08/wallHeatFlux'
    def zero(path,*a,**kw):
        text=actual_read(path,*a,**kw)
        return re.sub(r'value\s+uniform\s+[^;]+;','value uniform 0;',text) if path==target else text
    with patch.object(Path,'read_text',zero):
        assert all(r['heatFluxIntoWallWm2']==0 for r in read_wall_heat(CASE,'4e-08',32)['rows'])
    return dict(actualFrames=4,actualRowsBitwiseUnchanged=True,negativeAndZeroHeatPreserved=True,rejections=failures)

def region_checks():
    """模块二：从原表面独立核对映射、Heron面积、反向壁法向与封闭拓扑。"""
    c=json.loads((BASE/'region01/contract.json').read_text())
    s=json.loads((ROOT/'outputs/v194_validation/rocket_design03/surface.json').read_text())
    p=np.array(s['vertices']);f=np.array(s['facesOneBased'])-1
    data=np.load(BASE/'region01/boundary.npz');v=data['points'];tri=data['faces'];tags=data['patch']
    original_v=data['originalVertexIDs'];original_f=data['originalFaceIDs'];wall=tags==0
    assert sha(BASE/'region01/boundary.npz')==c['dataSHA256']
    for path,h in c['sourceSHA256'].items():assert sha(ROOT/path.replace('\\','/'))==h
    assert np.array_equal(v[original_v>0],p[original_v[original_v>0]-1])
    assert np.array_equal(original_f[wall],np.array(c['wallOriginalFaceIDs']))
    assert len(np.unique(original_f[wall]))==7838 and np.all(original_f[~wall]==0)
    assert np.array_equal(original_v[tri[wall]]-1,f[original_f[wall]-1][:,[0,2,1]])
    def heron(x):
        lengths=np.linalg.norm(x-np.roll(x,1,axis=1),axis=2);sp=lengths.sum(axis=1)/2
        return np.sqrt(sp*np.prod(sp[:,None]-lengths,axis=1))
    full=heron(p[f]);local=heron(v[tri[wall]])
    err=float(max(abs(full.sum()-c['fullOriginalAreaM2']),abs(local.sum()-c['wallAreaM2'])))
    assert err<1e-13 and len(f)==19752
    remainder=np.setdiff1d(np.arange(1,len(f)+1),original_f[wall])
    assert len(remainder)==c['remainingOriginalFaceCount']==11914
    assert abs(full[remainder-1].sum()-c['remainingOriginalAreaM2'])<1e-13
    directed=Counter((int(a),int(b)) for row in tri for a,b in zip(row,np.roll(row,-1)))
    assert all(n==1 and directed[(b,a)]==1 for (a,b),n in directed.items())
    assert len(v)-len(directed)//2+len(tri)==2
    # 边界是一个连通分量；所有点都由实际三角面使用。
    adjacency={i:set() for i in range(len(v))}
    for a,b in directed:adjacency[a].add(b)
    seen={0};todo=[0]
    while todo:
        for other in adjacency[todo.pop()]-seen:seen.add(other);todo.append(other)
    assert len(seen)==len(v)
    q=v[tri];sf=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])/2
    volume=float(np.einsum('ij,ij->i',q[:,0],sf).sum()/3)
    assert abs(volume/c['volumeM3']-1)<1e-14
    assert not c['executionAuthorized'] and not c['physicalQualification']
    assert not c['artificialPhysicalCap'] and c['newQualifiedHeatAreaM2']==0
    return dict(originalWallFaces=7838,remainingWallFaces=11914,originalWallAreaM2=float(local.sum()),
                independentHeronMaxAreaDifferenceM2=err,preparedAreaFraction=float(local.sum()/full.sum()),
                connectedClosedOrientedBoundary=True,volumeM3=volume,qualifiedHeatAreaM2=0)

def main():
    start=time.monotonic();out=BASE/'final_readonly01';out.mkdir(exist_ok=False)
    inputs=[p for case in ['thermal02','thermal03','region01'] for p in (BASE/case).rglob('*') if p.is_file()]
    before={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    shutil.copyfile(__file__,out/'executed_test.py');shutil.copyfile(ROOT/'scripts/readOpenFoamWallHeat.py',out/'executed_reader.py')
    try:
        reader=reader_checks();region=region_checks()
        protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
        assert len(protected)==29 and all(sha(ROOT/n.replace('\\','/'))==h for n,h in protected.items())
        assert all(sha(ROOT/n)==h for n,h in before.items())
        result=dict(passed=True,reader=reader,region=region,protectedFilesUnchanged=29,inputSHA256=before,
                    inputsUnchanged=True,CFDRuns=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
        assert result['elapsedSeconds']<=60
        save(out/'report.json',result)
        print(json.dumps(dict(passed=True,rejections=len(reader['rejections']),region=region,seconds=result['elapsedSeconds'])))
    except Exception as error:
        save(out/'report.json',dict(passed=False,error=str(error),elapsedSeconds=time.monotonic()-start));raise

if __name__=='__main__':main()
