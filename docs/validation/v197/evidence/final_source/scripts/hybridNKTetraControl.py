"""真四面体压缩—膨胀控制：成熟复数NK准备与唯一有限任务。"""
from pathlib import Path
import argparse,json,shutil,time
from collections import Counter
import numpy as np
from readEilmerNative import sha,read_grid,compare_source,read_field,need
from prepareHybridRocketEulerNative import local,save
from hybridRocketEuler import SavedFieldGuard

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v197_validation'
OLD=ROOT/'outputs/v191_validation/control01';CASE=BASE/'tetra01'
PACKAGE=ROOT/'tools/runtime/eilmer-nk-geometric-v1'
PROFILE='tetra-cea-compression-expansion-nk-120-v1'

NK='''-- 模块三：原生全表面谱半径、物理性检查和残差回溯；复数导数。
config.extrema_clipping=false
config.thermo_interpolator="rhop"
config.unstructured_limiter="hvenkat"
config.smooth_limiter_coeff=1.0
config.apply_unstructured_limiter_stagnation_point_filter=true
NewtonKrylovGlobalConfig{
 number_of_phases=1, use_preconditioner=true, preconditioner_perturbation=1.0e-50,
 preconditioner="ilu", ilu_fill=0, max_newton_steps=120, max_consecutive_bad_steps=10,
 stop_on_relative_residual=1.0e-8, frechet_derivative_perturbation=1.0e-50,
 use_scaling=true, max_linear_solver_iterations=50, max_linear_solver_restarts=0,
 inviscid_cfl_only=true, use_line_search=true, line_search_order=3,
 use_scaling_in_line_search=true, use_physicality_check=true,
 allowable_relative_mass_change=0.5, min_relaxation_factor_for_update=0.1,
 min_relaxation_factor_for_cfl_growth=0.5,
 number_of_steps_for_setting_reference_residuals=0, steps_between_status=1,
 write_loads=false, total_snapshots=8, steps_between_snapshots=20, steps_between_diagnostics=1
}
NewtonKrylovPhase:new{
 frozen_preconditioner=true, use_adaptive_preconditioner=true, steps_between_preconditioner_update=5,
 linear_solve_tolerance=0.01, residual_interpolation_order=2, jacobian_interpolation_order=1,
 use_residual_smoothing=true, use_auto_cfl=true, use_local_timestep=true,
 threshold_relative_residual_for_cfl_growth=0.99, start_cfl=0.1, max_cfl=100.0,
 auto_cfl_exponent=1.0, limit_on_cfl_decrease_ratio=1.0
}
'''

def design():
    """模块一：一致六分解，按实际四面体外边界继承标签，不用最近邻猜归属。"""
    CASE.mkdir(parents=True,exist_ok=False)
    lines=iter((OLD/'mesh.su2').read_text().splitlines())
    assert next(lines)=='NDIME= 3';nc=int(next(lines).split('=')[1])
    cells=[]
    for _ in range(nc):
        row=list(map(int,next(lines).split()));assert row[0]==12;cells.append(row[1:9])
    nv=int(next(lines).split('=')[1]);points=np.array([[float(x) for x in next(lines).split()[:3]] for _ in range(nv)])
    markers={}
    for _ in range(int(next(lines).split('=')[1])):
        tag=next(lines).split('=')[1].strip();n=int(next(lines).split('=')[1]);markers[tag]=[]
        for _ in range(n):
            row=list(map(int,next(lines).split()));assert row[0]==9;markers[tag].append(row[1:5])
    assert (nc,nv)==(2400,3843)
    split=np.array([[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]])
    tetra=np.array(cells)[:,split].reshape(-1,4)
    xyz=points[tetra];det=np.linalg.det((xyz[:,1:]-xyz[:,:1]).transpose(0,2,1))
    neg=det<0;tetra[neg,0],tetra[neg,1]=tetra[neg,1].copy(),tetra[neg,0].copy()
    assert np.all(abs(det)>0)
    actual_volume=float(abs(det).sum()/6)
    analytical=.04*(.6-.24*np.tan(np.deg2rad(10)))
    assert abs(actual_volume/analytical-1)<2e-13
    faces={};count=Counter();signs={}
    from itertools import combinations
    source_tri={}
    for tag,quads in markers.items():
        for quad in quads:
            for tri in combinations(quad,3):
                key=tuple(sorted(tri));assert key not in source_tri;source_tri[key]=tag
    for cell in tetra:
        for order in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]:
            tri=tuple(cell[order]);key=tuple(sorted(tri));count[key]+=1;faces[key]=tri
            inversions=sum(tri[i]>tri[j] for i in range(3) for j in range(i+1,3))
            signs[key]=signs.get(key,0)+(-1 if inversions%2 else 1)
    assert all(n in (1,2) for n in count.values())
    assert all(signs[key]==0 for key,n in count.items() if n==2)
    boundaries={tag:[] for tag in markers}
    for key,n in count.items():
        if n==1:
            assert key in source_tri;boundaries[source_tri[key]].append(faces[key])
    assert all(len(boundaries[k])==2*len(markers[k]) for k in markers)
    with (CASE/'mesh.su2').open('w') as f:
        f.write('NDIME= 3\nNELEM= %d\n'%len(tetra))
        for i,row in enumerate(tetra):f.write('10 '+' '.join(map(str,row))+' '+str(i)+'\n')
        f.write('NPOIN= %d\n'%nv)
        for i,row in enumerate(points):f.write(' '.join(format(x,'.17g') for x in row)+' '+str(i)+'\n')
        f.write('NMARK= %d\n'%len(boundaries))
        for tag,rows in boundaries.items():
            f.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n'%(tag,len(rows)))
            for row in rows:f.write('5 '+' '.join(map(str,row))+'\n')
    for name in ['grid.lua','rah-air.gas']:shutil.copyfile(OLD/name,CASE/name)
    job=(OLD/'job.lua').read_text(encoding='utf-8').split('config.gasdynamic_update_scheme=')[0]
    assert job.count('config.solver_mode = "time-marching"')==1
    job=job.replace('config.solver_mode = "time-marching"','config.solver_mode = "steady"')
    (CASE/'job.lua').write_text(job+NK,encoding='utf-8')
    wall=np.array(boundaries['wall']);tri=points[wall]
    area=float(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()/2)
    save(CASE/'surface.json',dict(vertices=points.tolist(),facesOneBased=(wall+1).tolist(),
        faceIDs=list(range(1,len(wall)+1)),identity='new tetrahedral control; not original rocket wall'))
    m=dict(schema='RAH.HybridNKTetraControl.v1',profile=PROFILE,cells=len(tetra),vertices=nv,faces=len(count),
        ghosts=sum(len(x) for x in boundaries.values()),volumeM3=actual_volume,wallTriangles=len(wall),wallAreaM2=area,
        runtimeSHA256=sha(PACKAGE/'runtime.json'),sourceHexMeshSHA256=sha(OLD/'mesh.su2'),
        filesSHA256={n:sha(CASE/n) for n in ['mesh.su2','grid.lua','rah-air.gas','job.lua','surface.json']},
        physicalQualification=False,executed=False)
    save(CASE/'manifest.json',m);return m

def prepare():
    """模块二：原生准备和全体配对/物性审核共用90秒；无求解推进。"""
    from runEilmerProcess import run
    start=time.monotonic();m=design();out=BASE/'preparation01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    commands=[]
    for i,cmd in enumerate([['lmr','prep-grids','--job=grid.lua'],['lmr','prep-sim','--job=job.lua']]):
        p=run(cmd,CASE,out/('process'+str(i)),90-(time.monotonic()-start),96,[Path(__file__),*CASE.glob('*.*')],4096,4096)
        commands.append(p);assert p['passed']
    grid=read_grid(CASE/'lmrsim/grid/grid-0000');compare_source(grid,CASE/'mesh.su2')
    assert len(grid['cells'])==14400 and np.all(grid['kinds']==10)
    f=read_field(CASE,0,14400);assert np.all(f['vol']>0) and abs(f['vol'].sum()/m['volumeM3']-1)<2e-13
    cfg=json.loads((CASE/'lmrsim/config').read_text());nk=json.loads((CASE/'lmrsim/nk-config').read_text())
    assert cfg['solver_mode']=='steady' and cfg['viscous'] is False and cfg['turbulence_model']=='none'
    for key in ['inviscid_cfl_only','use_line_search','use_physicality_check']:assert nk[key] is True
    assert nk['frechet_derivative_perturbation']==1e-50 and nk['max_newton_steps']==120
    result=dict(passed=True,case=str(CASE),geometry=m,commands=commands,
        nativeInputSHA256={str(p):sha(p) for p in CASE.rglob('*') if p.is_file()},
        elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    save(out/'report.json',result);print(json.dumps(dict(passed=True,cells=m['cells'],elapsedSeconds=result['elapsedSeconds'])))
    return True

def execute():
    """模块四：一次原生复数NK控制；完整场守卫独立于原生物理性检查。"""
    from runEilmerProcess import run
    start=time.monotonic();p=BASE/'preparation01/report.json';r=json.loads(p.read_text())
    assert r['passed'] and r['geometry']['profile']==PROFILE
    runtime=PACKAGE/'runtime.json';pkg=json.loads(runtime.read_text());assert sha(runtime)==r['geometry']['runtimeSHA256']
    assert pkg['profile']=='rah-nk-geometric-v1' and pkg['numberType']=='complex'
    for path,h in {**pkg['filesSHA256'],**r['nativeInputSHA256']}.items():assert sha(local(path))==h
    with (CASE/'nk-attempt.json').open('x') as f:json.dump(dict(profile=PROFILE,automaticContinuation=0),f)
    out=BASE/'run01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    guard=SavedFieldGuard(CASE,14400,out/'saved_field_guard.json')
    locks=[Path(__file__),ROOT/'scripts/hybridRocketEuler.py',p,runtime,ROOT/'docs/validation/v197/ENTRY_DECISION.md',
           *pkg['filesSHA256'],*r['nativeInputSHA256']]
    process=run([pkg['binary'],'--max-cpus=2','--max-wall-clock=240'],CASE,out/'native',
        260-(time.monotonic()-start),96,locks,4096,4096,checkpoint_guard=guard)
    result=dict(profile=PROFILE,processPassed=process['passed'],process=process,elapsedSeconds=time.monotonic()-start,
        outputSHA256={str(p):sha(p) for p in (CASE/'lmrsim').rglob('*') if p.is_file()},physicalQualification=False)
    save(out/'report.json',result);print(json.dumps(dict(passed=process['passed'],exit=process['nativeExitCode'],
        stop=process['stopReason'],elapsedSeconds=result['elapsedSeconds'])));return process['passed']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','execute']);a=p.parse_args()
    raise SystemExit(0 if (prepare() if a.stage=='prepare' else execute()) else 1)
