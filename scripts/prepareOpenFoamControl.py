"""复用准确SU2四面体连接生成原生polyMesh及NASA7无黏控制；不启动CFD。"""
from pathlib import Path
import argparse,json,time
import numpy as np
from openFoamAirThermo import mixture_dictionary,contract,evaluate
from hashlib import sha256

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v198_validation'
CASE=BASE/'control01';OLD=ROOT/'outputs/v197_validation/tetra01'
def sha(p):return sha256(Path(p).read_bytes()).hexdigest()
def save(p,r):Path(p).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def header(name,kind='dictionary'):
    return 'FoamFile { format ascii; class '+kind+'; object '+name+'; }\n'


def read_tetra(path):
    """模块一：仅接受明确3D四面体/三角面及顺序编号，不隐式截断其他网格。"""
    lines=iter(Path(path).read_text().splitlines());assert next(lines)=='NDIME= 3'
    n=int(next(lines).split('=')[1]);cells=[]
    for i in range(n):
        row=list(map(int,next(lines).split()));assert len(row)==6 and row[0]==10 and row[-1]==i
        cells.append(row[1:5])
    nv=int(next(lines).split('=')[1]);pts=[]
    for i in range(nv):
        row=next(lines).split();assert len(row)==4 and int(row[-1])==i;pts.append(list(map(float,row[:3])))
    markers={}
    for _ in range(int(next(lines).split('=')[1])):
        tag=next(lines).split('=')[1].strip();assert tag not in markers
        count=int(next(lines).split('=')[1]);markers[tag]=[]
        for i in range(count):
            row=list(map(int,next(lines).split()));assert len(row)==4 and row[0]==5;markers[tag].append(row[1:])
    assert next(lines,None) is None
    p,c=np.array(pts),np.array(cells)
    assert np.isfinite(p).all() and c.min()>=0 and c.max()<nv
    assert len({tuple(sorted(x)) for x in cells})==len(cells)
    return p,c,markers


def topology(p,c,markers):
    """模块二：按连接精确配对，面法向指向owner外侧；边界保持原标签与行序。"""
    centers=p[c].mean(axis=1);faces={}
    volumes=np.abs(np.linalg.det((p[c][:,1:]-p[c][:,:1]).transpose(0,2,1)))/6
    assert np.all(volumes>0)
    for owner,cell in enumerate(c):
        for order in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]:
            tri=cell[order].copy();key=tuple(sorted(tri));xyz=p[tri]
            if np.dot(np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0]),xyz.mean(axis=0)-centers[owner])<0:tri=tri[::-1]
            if key not in faces:faces[key]=[tri,owner,None]
            else:
                assert faces[key][2] is None,'非流形面';faces[key][2]=owner
    internal=sorted([v for v in faces.values() if v[2] is not None],key=lambda v:(v[1],v[2]))
    boundary=[];patches=[];seen=set()
    for tag,rows in markers.items():
        start=len(internal)+len(boundary)
        for original,row in enumerate(rows):
            key=tuple(sorted(row));assert key not in seen;seen.add(key)
            assert key in faces and faces[key][2] is None,'边界不属于实际外表面'
            boundary.append(faces[key])
        patches.append(dict(name=tag,startFace=start,nFaces=len(rows)))
    assert seen=={key for key,v in faces.items() if v[2] is None},'边界遗漏'
    ordered=internal+boundary;tri=np.array([v[0] for v in ordered]);owner=np.array([v[1] for v in ordered])
    neighbor=np.array([v[2] for v in internal]);xyz=p[tri]
    sf=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])/2;fc=xyz.mean(axis=1)
    assert np.all(np.einsum('ij,ij->i',sf,fc-centers[owner])>0)
    assert np.all(owner[:len(neighbor)]<neighbor)
    assert np.all(np.einsum('ij,ij->i',sf[:len(neighbor)],centers[neighbor]-centers[owner[:len(neighbor)]])>0)
    closure=np.zeros((len(c),3));np.add.at(closure,owner,sf);np.add.at(closure,neighbor,-sf[:len(neighbor)])
    assert np.max(np.abs(closure))<1e-14
    return tri,owner,neighbor,patches,dict(cells=len(c),points=len(p),faces=len(tri),internalFaces=len(neighbor),
        volumeM3=float(volumes.sum()),minimumCellVolumeM3=float(volumes.min()),maxCellAreaClosureM2=float(abs(closure).max()))


def prepare(repair=False):
    """模块三：只在新目录写入完整原生配置与原面映射，默认endTime=0。"""
    start=time.monotonic();CASE=BASE/('control03' if repair else 'control01');CASE.mkdir(parents=True,exist_ok=False)
    m=json.loads((OLD/'manifest.json').read_text());assert sha(OLD/'mesh.su2')==m['filesSHA256']['mesh.su2']
    p,c,markers=read_tetra(OLD/'mesh.su2');tri,owner,neighbor,patches,geo=topology(p,c,markers)
    assert geo['cells']==14400 and geo['points']==3843 and geo['faces']==31520 and len(markers['wall'])==240
    repair_record=None
    if repair:
        # 模块四：只细分多边界角单元；重心到四个原面，外三角连接不变。
        internal_counts=np.bincount(np.r_[owner[:len(neighbor)],neighbor],minlength=len(c))
        selected=np.flatnonzero(internal_counts<3);assert len(selected)==164
        newpts=p[c[selected]].mean(axis=1);newcells=[];parents=[];indices=set(selected.tolist());idx=0
        for i,cell in enumerate(c):
            if i not in indices:newcells.append(cell.tolist());parents.append(i);continue
            center=len(p)+idx;idx+=1
            for face in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]:
                newcells.append([center,*cell[face].tolist()]);parents.append(i)
        p=np.vstack([p,newpts]);c=np.array(newcells)
        oldvol=geo['volumeM3'];tri,owner,neighbor,patches,geo=topology(p,c,markers)
        assert geo['cells']==14892 and len(p)==4007 and abs(geo['volumeM3']/oldvol-1)<2e-13
        repair_record=dict(selectedOriginalCellsZeroBased=selected.tolist(),childParentCellZeroBased=parents,
            method='one centroid four-split only where internal-face count <3',originalBoundaryExact=True,
            originalCheckMeshOverallPassed=False)
        save(CASE/'repair.json',repair_record)
    mesh=CASE/'constant/polyMesh';mesh.mkdir(parents=True);(CASE/'system').mkdir();(CASE/'0').mkdir()
    def listing(name,kind,rows):
        (mesh/name).write_text(header(name,kind)+str(len(rows))+'\n(\n'+'\n'.join(rows)+'\n)\n')
    listing('points','vectorField',['('+' '.join(format(x,'.17g') for x in row)+')' for row in p])
    listing('faces','faceList',['3('+' '.join(map(str,row))+')' for row in tri])
    listing('owner','labelList',list(map(str,owner)));listing('neighbour','labelList',list(map(str,neighbor)))
    listing('boundary','polyBoundaryMesh',[q['name']+' { type '+('wall' if q['name'] in ['wall','span'] else 'patch')+
        '; nFaces '+str(q['nFaces'])+'; startFace '+str(q['startFace'])+'; }' for q in patches])
    thermo=header('physicalProperties')+'thermoType { type hePsiThermo; mixture pureMixture; transport const; thermo janaf; equationOfState perfectGas; specie specie; energy sensibleInternalEnergy; }\nmixture {\n'+mixture_dictionary()+'}\n'
    (CASE/'constant/physicalProperties').write_text(thermo)
    g=evaluate(5000,300);u=-3*float(g['soundSpeedMps']);values={'p':'5000','T':'300','U':'(%.17g 0 0)'%u}
    dims={'p':'1 -1 -2 0 0 0 0','T':'0 0 0 1 0 0 0','U':'0 1 -1 0 0 0 0'}
    for name,value in values.items():
        rows=[]
        for patch in patches:
            tag=patch['name']
            if tag=='inflow':bc='type fixedValue; value uniform '+value+';'
            elif tag=='outer':bc='type freestream; freestreamValue uniform '+value+'; value uniform '+value+';'
            elif tag in ['wall','span'] and name=='U':bc='type slip;'
            else:bc='type zeroGradient;'
            rows.append(tag+' { '+bc+' }')
        (CASE/'0'/name).write_text(header(name,'volVectorField' if name=='U' else 'volScalarField')+
            'dimensions ['+dims[name]+'];\ninternalField uniform '+value+';\nboundaryField {\n'+'\n'.join(rows)+'\n}\n')
    (CASE/'system/controlDict').write_text(header('controlDict')+'''solver shockFluid;
startFrom startTime; startTime 0; stopAt endTime; endTime 0; deltaT 1;
writeControl timeStep; writeInterval 20; purgeWrite 0; writeFormat ascii;
writePrecision 17; timeFormat general; timePrecision 12; writeCompression off;
runTimeModifiable false; adjustTimeStep no;
functions { energy { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e); writeControl writeTime; } }
''')
    (CASE/'system/fvSchemes').write_text(header('fvSchemes')+'''fluxScheme Kurganov;
ddtSchemes { default localEuler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; reconstruct(rho) vanLeer; reconstruct(U) vanLeerV; reconstruct(T) vanLeer; }
snGradSchemes { default corrected; }
''')
    (CASE/'system/fvSolution').write_text(header('fvSolution')+'''solvers {
 "rho.*" { solver diagonal; }
 "U.*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol 0; }
 "e.*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0; }
}
PIMPLE { maxCo 0.8; rDeltaTSmoothingCoeff 0.02; }
''')
    np.savez(CASE/'exact_mesh.npz',points=p,cells=c,faces=tri,owner=owner,neighbor=neighbor)
    save(CASE/'patches.json',patches);save(CASE/'thermo_contract.json',contract())
    wall=next(x for x in patches if x['name']=='wall')
    ids=np.arange(wall['startFace'],wall['startFace']+wall['nFaces'])
    assert np.array_equal(np.sort(tri[ids],axis=1),np.sort(np.array(markers['wall']),axis=1))
    save(CASE/'wall_identity.json',dict(originalFaceIDs=list(range(1,len(ids)+1)),foamFaceIDsZeroBased=ids.tolist(),
        originalConnections=markers['wall'],identity='240 control triangles, NOT original 19752 rocket faces',exact=True))
    report=dict(passed=True,geometry=geo,sourceMeshSHA256=sha(OLD/'mesh.su2'),
        patches=patches,condition=dict(pPa=5000,TK=300,Mach=3,velocityMps=[u,0,0]),
        inputSHA256={str(x):sha(x) for x in CASE.rglob('*') if x.is_file()},
        elapsedSeconds=time.monotonic()-start,physicalQualification=False,CFDRuns=0)
    save(BASE/('preparation_repaired.json' if repair else 'preparation.json'),report);print(json.dumps(dict(passed=True,geometry=geo,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--repair-corner-stencil',action='store_true');args=parser.parse_args()
    prepare(args.repair_corner_stencil)
