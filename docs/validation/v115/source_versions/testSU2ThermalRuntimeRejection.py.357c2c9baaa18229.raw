"""实际调用受控引擎验证遗漏标记/文件和覆盖不足均失败，不执行科学长任务。"""
from pathlib import Path
import hashlib,json,os,subprocess,time,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from su2ThermalRuntime import resolve
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=ROOT/'outputs/v115_validation';base=out/'actual01/native_constant'
    dest=out/'native_rejections';assert not dest.exists();dest.mkdir()
    runtime=ROOT/'tools/runtime/su2/8.6.0-rah-axial-wall-v1/runtime.json'
    exe=Path(json.loads(runtime.read_text())['executable']);config=(base/'case.cfg').read_text(encoding='utf-8')
    reports=[]
    for name,change,table,reason in [
        ('missing_marker',lambda c:c.replace('MARKER_PYTHON_CUSTOM= (wall)',''),(base/'wall_temperature.csv').read_bytes(),'exactly one custom wall marker'),
        ('missing_table',lambda c:c,None,'cannot be read'),
        ('incomplete_coverage',lambda c:c,b'RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K\n0.1,100\n0.2,110\n','does not cover wall node X')]:
        folder=dest/name;folder.mkdir();cfg=folder/'case.cfg';cfg.write_text(change(config),encoding='utf-8')
        os.link(str(base/'mesh.su2'),str(folder/'mesh.su2'))
        profile=folder/'wall_temperature.csv'
        if table is not None:profile.write_bytes(table)
        lock={str(p):sha(p) for p in [Path(__file__).resolve(),exe,runtime,cfg,folder/'mesh.su2']+([profile] if table else [])}
        command=[str(exe),str(cfg),'--threads','2'];start=time.monotonic()
        with (folder/'engine_stdout.log').open('wb') as f:
            proc=subprocess.Popen(command,cwd=str(folder),stdout=f,stderr=f,creationflags=subprocess.CREATE_NO_WINDOW,
                env=dict(os.environ,RAH_WALL_TEMPERATURE_FILE=str(profile),OMP_NUM_THREADS='2'))
            timed=False
            try:code=proc.wait(timeout=20)
            except subprocess.TimeoutExpired:timed=True;proc.kill();code=proc.wait()
        text=(folder/'engine_stdout.log').read_text(encoding='utf-8',errors='replace')
        row=dict(name=name,status='terminal',processID=proc.pid,exitCode=code,timedOut=timed,
            elapsedSeconds=time.monotonic()-start,command=command,inputSHA256=lock,
            inputsStable=all(sha(Path(p))==h for p,h in lock.items()),expectedMessage=reason,
            expectedMessageFound=reason in text,solverIterationStarted='Begin Solver' in text)
        (folder/'process.json').write_text(json.dumps(row,indent=2)+'\n');reports.append(row)
        assert code!=0 and not timed and row['inputsStable'] and reason in text
        assert not (folder/'surface_flow.vtk').exists()
    # 监督器也须在启动前拒绝输入清单与配置失配。
    rejected=[]
    for name,c in [('missing_marker',config.replace('MARKER_PYTHON_CUSTOM= (wall)','')),
                   ('multigrid',config.replace('MGLEVEL= 0','MGLEVEL= 1')),
                   ('unsteady',config+'\nTIME_DOMAIN= YES\n')]:
        try:resolve(base,c)
        except AssertionError:rejected.append(name)
    assert len(rejected)==3
    (out/'runtime_rejections.json').write_text(json.dumps(dict(native=reports,supervisor=rejected,physicsQualified=False),indent=2)+'\n')
    print('Three actual native refusals and three prelaunch configuration refusals passed.')
if __name__=='__main__':main()
