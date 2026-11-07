"""真实原19752面OpenFOAM入口：精确转换、全壁映射及相容NASA7初场。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save,header,read_tetra,topology
from openFoamAirThermo import evaluate
from hybridEdgeSampling import map_original_wall

BASE=ROOT/'outputs/v200_validation';SOURCE=ROOT/'outputs/v194_validation/rocket_design03'


def main():
    start=time.monotonic();out=BASE/'rocket_preparation02';out.mkdir(exist_ok=False);case=BASE/'rocket02';case.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');shutil.copyfile(ROOT/'scripts/prepareOpenFoamControl.py',out/'topology_source.py')
    assert json.loads((BASE/'restart_audit02/report.json').read_text())['passed']
    m=json.loads((SOURCE/'manifest.json').read_text());surface=json.loads((SOURCE/'surface.json').read_text())
    for name in ['mesh.su2','surface.json']:assert sha(SOURCE/name)==m['filesSHA256'][name]
    p,c,markers=read_tetra(SOURCE/'mesh.su2',implicit_element_ids=True);tri,owner,neighbor,patches,geo=topology(p,c,markers)
    assert len(c)==64982 and len(p)==18027 and set(markers)=={'wall','farfield'}
    assert len(markers['wall'])==19752 and len(markers['farfield'])==6572
    # 模块一：必要时仅修复原生内部面方向缺秩单元；全部原边界及旧网格保留。
    counts=np.bincount(np.r_[owner[:len(neighbor)],neighbor],minlength=len(c));selected=np.flatnonzero(counts<3)
    parents=list(range(len(c)))
    if len(selected):
        oldp=p.copy();newcells=[];parents=[];lookup={int(k):len(p)+j for j,k in enumerate(selected)}
        for i,cell in enumerate(c):
            if i not in lookup:newcells.append(cell);parents.append(i);continue
            for face in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]:newcells.append([lookup[i],*cell[face]]);parents.append(i)
        p=np.vstack([p,p[c[selected]].mean(axis=1)]);c=np.array(newcells)
        tri,owner,neighbor,patches,newgeo=topology(p,c,markers)
        assert abs(newgeo['volumeM3']/geo['volumeM3']-1)<2e-13;geo=newgeo
    mesh=case/'constant/polyMesh';mesh.mkdir(parents=True);(case/'system').mkdir();(case/'0').mkdir()
    def listing(name,kind,rows):
        (mesh/name).write_text(header(name,kind)+str(len(rows))+'\n(\n'+'\n'.join(rows)+'\n)\n')
    listing('points','vectorField',['('+' '.join(format(x,'.17g') for x in row)+')' for row in p])
    listing('faces','faceList',['3('+' '.join(map(str,row))+')' for row in tri])
    listing('owner','labelList',list(map(str,owner)));listing('neighbour','labelList',list(map(str,neighbor)))
    listing('boundary','polyBoundaryMesh',[q['name']+' { type '+('wall' if q['name']=='wall' else 'patch')+'; nFaces '+str(q['nFaces'])+'; startFace '+str(q['startFace'])+'; }' for q in patches])
    np.savez(case/'exact_mesh.npz',points=p,cells=c,faces=tri,owner=owner,neighbor=neighbor)
    save(case/'patches.json',patches)
    # 模块二：复用已有精确连接双射；此处采样位置为四面体体积质心，尚无原生发展场。
    wall=next(q for q in patches if q['name']=='wall');ids=np.arange(wall['startFace'],wall['startFace']+wall['nFaces'])
    grid=dict(points=p,cells=c,faces=tri,boundaries={'wall':dict(ids=ids,owners=owner[ids],signs=np.ones(len(ids),int))})
    mapping=map_original_wall(surface,grid,dict(originalNodes=markers['wall'],nativeRows=list(range(len(ids)))),p[c].mean(axis=1))
    mapping['positionKind']='exact-tetra-volume-centroid; native field export pending'
    xyz=p[tri[ids]];area=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)/2
    assert abs(area.sum()-m['geometry']['wallAreaM2'])<2e-13
    mapping.update(originalAreaM2=float(area.sum()),areaM2=area.tolist(),originalConnections=surface['facesOneBased'])
    save(case/'wall_mapping.json',mapping)
    # 模块三：固定Mach/方向，按本NASA7声速重算速度；不挪用CEA速度当同物性。
    condition=m['condition'];gas=evaluate(condition['pressurePa'],condition['temperatureK'])
    u=condition['Mach']*float(gas['soundSpeedMps'])*np.array(condition['velocityDirection'])
    old=json.loads((ROOT/'outputs/v194_validation/rocket01/native-input-state.json').read_text())
    values={'p':format(condition['pressurePa'],'.17g'),'T':format(condition['temperatureK'],'.17g'),'U':'('+' '.join(format(x,'.17g') for x in u)+')'}
    dims={'p':'1 -1 -2 0 0 0 0','T':'0 0 0 1 0 0 0','U':'0 1 -1 0 0 0 0'}
    for name,value in values.items():
        wallbc='type slip;' if name=='U' else 'type zeroGradient;'
        (case/'0'/name).write_text(header(name,'volVectorField' if name=='U' else 'volScalarField')+
            'dimensions ['+dims[name]+'];\ninternalField uniform '+value+';\nboundaryField {\nwall { '+wallbc+' }\nfarfield { type freestream; freestreamValue uniform '+value+'; value uniform '+value+'; }\n}\n')
    template=BASE/'restart01'
    for name in ['physicalProperties']:shutil.copyfile(template/'constant'/name,case/'constant'/name)
    for name in ['fvSchemes','fvSolution']:shutil.copyfile(template/'system'/name,case/'system'/name)
    shutil.copyfile(template/'thermo_contract.json',case/'thermo_contract.json')
    control=(template/'system/controlDict').read_text().replace('startTime 0.002;','startTime 0;').replace('endTime 0.003;','endTime 0;')
    (case/'system/controlDict').write_text(control)
    report=dict(passed=True,profile='of13-nasa7-original-wall-rocket-entry-v1',geometry=geo,
        originalWallFaces=19752,originalWallAreaM2=float(area.sum()),originalFarfieldFaces=6572,
        repairedOriginalCellsZeroBased=selected.tolist(),newCellCount=len(c),originalWallExact=True,
        minimumWallNormalDistanceM=min(mapping['normalDistanceM']),maximumWallNormalDistanceM=max(mapping['normalDistanceM']),
        condition=condition,gasSoundSpeedMps=float(gas['soundSpeedMps']),velocityMps=u.tolist(),
        oldCEAVelocityMps=[old[k] for k in ['velx','vely','velz']],oldCEASoundSpeedMps=old['a'],
        speedRelativeChangeFromCEA=float(gas['soundSpeedMps']/old['a']-1),
        boundaryMethod='stock freestream: fixed incoming state / zeroGradient outgoing; slip wall',
        endTime=0,executionAuthorized=False,physicalQualification=False,newQualifiedHeatAreaM2=0,
        temperatureIdentity='initial fluid temperature, not solved solid Tw or qualified edge',
        sourceSHA256={str(q):sha(q) for q in [SOURCE/'manifest.json',SOURCE/'mesh.su2',SOURCE/'surface.json',BASE/'restart_audit02/report.json',ROOT/'docs/validation/v200/ROCKET_ENTRY_PLAN.md']},
        inputSHA256={str(q):sha(q) for q in case.rglob('*') if q.is_file()},elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','geometry','repairedOriginalCellsZeroBased','velocityMps','speedRelativeChangeFromCEA','elapsedSeconds']}))


if __name__=='__main__':main()
