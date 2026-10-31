"""工作台变物性后台流程：初态转换、原生求解、EOS审核分别记账，不伪造终态。"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, json, os, subprocess, sys, time
from runEmbeddedSU2 import ROOT, sha, atomic_write_json


class PipelineStopped(RuntimeError):
    pass


class Pipeline:
    def __init__(self, request_path):
        self.request_path = Path(request_path).resolve()
        self.request = json.loads(self.request_path.read_text(encoding='utf-8-sig'))
        self.folder = Path(self.request['workDirectory']).resolve()
        self.folder.relative_to(ROOT / 'outputs')
        assert self.request_path.parent == self.folder and self.request['jobToken']
        assert ('propertyInitialization' in self.request)!=('turbulenceInitialization' in self.request)
        self.sst='turbulenceInitialization' in self.request
        self.spec = self.request['turbulenceInitialization' if self.sst else 'propertyInitialization']
        if self.sst:assert self.spec['model']=='SST-2003m'
        else:assert self.spec['model'] == 'HEOS-Air' and self.spec['runtimeFlavor'] == 'rah-coolprop-recon-v3'
        self.out = self.folder / 'workflow'
        self.out.mkdir(exist_ok=False)
        self.file = self.out / 'pipeline_process.json'
        self.start = time.monotonic()
        paths = [self.request_path, Path(__file__).resolve()]
        paths += [ROOT / n for n in ('scripts/prepareSU2CoolPropCase.py', 'scripts/auditSU2CoolPropState.py',
                                    'scripts/runEmbeddedSU2.py', 'scripts/su2ThermalRuntime.py')]
        if self.sst:
            paths += [ROOT/n for n in ('scripts/prepareSU2SSTCase.py','scripts/su2SSTInterface.py','scripts/auditSU2SSTState.py')]
        self.locks = {str(p): sha(p) for p in paths}
        self.state = dict(status='starting', stage='preparation', supervisorProcessID=os.getpid(),
                          jobToken=self.request['jobToken'], startedUTC=datetime.now(timezone.utc).isoformat(),
                          phases=[], engineStarted=False, cancelled=False, timedOut=False,
                          inputSHA256=self.locks, physicsQualified=False, newQualifiedArea=0)
        self.stage_start = self.start
        self.deadline = self.start

    def save(self):
        self.state.update(elapsedSeconds=time.monotonic()-self.start,
                          phaseElapsedSeconds=time.monotonic()-self.stage_start)
        atomic_write_json(self.file, self.state)

    def cancelled(self):
        try:
            return json.loads((self.folder/'engine_cancel.json').read_text(encoding='utf-8-sig')).get('jobToken') == self.request['jobToken']
        except (ValueError, OSError):
            return False

    def check(self):
        self.state['cancelled'] = self.cancelled()
        self.state['timedOut'] = time.monotonic() >= self.deadline
        if self.state['cancelled'] or self.state['timedOut']:
            raise PipelineStopped('工作台阶段取消' if self.state['cancelled'] else '工作台阶段超过独立预算')

    def run_property(self, out, name, command, cwd, env, seconds):
        """物性程序为直接子进程；只结束本次拥有的进程，保存每项真实退出码。"""
        self.check()
        receipt=out/(name+'_process.json');assert not receipt.exists()
        start=time.monotonic();limit=min(float(seconds),self.deadline-start)
        state=dict(status='starting',command=list(map(str,command)),cwd=str(cwd),maximumSeconds=limit)
        with (out/(name+'.log')).open('wb') as log:
            proc=subprocess.Popen(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,
                                  creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            state.update(status='running',processID=proc.pid)
            atomic_write_json(receipt,state)
            try:
                while proc.poll() is None:
                    self.check()
                    if time.monotonic()-start >= limit:
                        self.state['timedOut']=True;raise PipelineStopped('物性查询超时')
                    self.save();time.sleep(.1)
            finally:
                if proc.poll() is None:proc.kill()
                code=proc.wait()
                state.update(status='terminal',exitCode=code,cancelled=self.state['cancelled'],
                             timedOut=self.state['timedOut'],elapsedSeconds=time.monotonic()-start)
                atomic_write_json(receipt,state)
        assert code==0, '物性程序未成功：'+str(receipt)

    def run_engine(self):
        """原生监督器继续拥有SU2；取消先交给它写真实终态，异常才回收此进程树。"""
        self.check()
        command=[sys.executable,'-X','utf8',str(ROOT/'scripts/runEmbeddedSU2.py'),'--request',str(self.request_path)]
        with (self.out/'native_supervisor.log').open('wb') as log:
            proc=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,
                                  creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            self.state.update(engineStarted=True,nativeSupervisorProcessID=proc.pid);self.save()
            force_deadline=None
            try:
                while proc.poll() is None:
                    if self.cancelled() or time.monotonic()>=self.deadline:
                        self.state['cancelled']=self.cancelled()
                        self.state['timedOut']=time.monotonic()>=self.deadline
                        atomic_write_json(self.folder/'engine_cancel.json',{'jobToken':self.request['jobToken']})
                        if force_deadline is None:force_deadline=time.monotonic()+10
                        if time.monotonic()>=force_deadline:
                            raise PipelineStopped('原生监督器未在取消宽限内退出')
                    self.save();time.sleep(.2)
            finally:
                if proc.poll() is None:
                    # PID来自刚刚创建的活进程；不按名称终止其他用户计算。
                    subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],stdout=log,stderr=log,
                                   creationflags=subprocess.CREATE_NO_WINDOW,check=False)
                code=proc.wait();self.state['nativeSupervisorExitCode']=code
            native=self.folder/'engine_process.json'
            assert native.exists(), '原生监督器没有产生引擎记录，不能伪造SU2终态'
            r=json.loads(native.read_text(encoding='utf-8-sig'))
            self.state.update(nativeReceiptSHA256=sha(native),cancelled=bool(r.get('cancelled')),
                              timedOut=bool(r.get('timedOut')))
            assert code==0 and r['status']=='terminal' and r['exitCode']==0 and not r['cancelled'] and not r['timedOut'], '原生计算未成功'
            assert all(r['inputsStable'].values()), '原生计算输入改变'

    def phase(self, name, seconds, action):
        assert 1 <= seconds <= 750
        self.stage_start=time.monotonic();self.deadline=self.stage_start+seconds
        self.state.update(status='running',stage=name,phaseBudgetSeconds=seconds);self.save()
        try:
            self.check();action();self.check()
        finally:
            self.state['phases'].append(dict(stage=name,maximumSeconds=seconds,
                                             elapsedSeconds=time.monotonic()-self.stage_start))

    def execute(self):
        from prepareSU2CoolPropCase import prepare
        from auditSU2CoolPropState import audit
        code=1
        try:
            if self.sst:
                command=[sys.executable,'-X','utf8',str(ROOT/'scripts/prepareSU2SSTCase.py'),
                         '--folder',str(self.folder),'--source',self.spec['sourceDirectory'],
                         '--intensity',str(self.spec['intensity']),'--viscosity-ratio',str(self.spec['viscosityRatio']),
                         '--prt',str(self.spec['Prt']),'--provenance',self.spec['provenance']]
                self.phase('preparation',180,lambda:self.run_property(self.out,'sst_preparation',command,ROOT,dict(os.environ),180))
            else:
                self.phase('preparation',180,lambda:prepare(self.folder,self.spec['sourceDirectory'],
                                                           self.spec['runtimeFlavor'],self.run_property))
            # 原生预算由旧监督器严格执行；额外30秒仅供启动/清场，不增加SU2预算。
            self.phase('solver',float(self.request['maximumSeconds'])+30,self.run_engine)
            if self.sst:
                target=self.folder/'sst_state_audit.json'
                command=[sys.executable,'-X','utf8',str(ROOT/'scripts/auditSU2SSTState.py'),'--folder',str(self.folder),'--output',str(target)]
                self.phase('eosAudit',120,lambda:self.run_property(self.out,'sst_audit',command,ROOT,dict(os.environ),120))
            else:
                self.phase('eosAudit',120,lambda:audit(self.folder,self.run_property));target=self.folder/'thermodynamics/audit.json'
            self.state['auditSHA256']=sha(target)
            code=0
        except Exception as ex:
            self.state['error']=type(ex).__name__+': '+str(ex)
        finally:
            stable={p:sha(Path(p))==h for p,h in self.locks.items()}
            if not all(stable.values()):code=2
            if self.state['cancelled']:code=125
            elif self.state['timedOut']:code=124
            self.state.update(status='terminal',exitCode=code,inputsStable=stable,
                              endedUTC=datetime.now(timezone.utc).isoformat())
            self.save()
        return code


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--request',required=True);args=parser.parse_args()
    pipeline=Pipeline(args.request);code=pipeline.execute()
    print(json.dumps({k:pipeline.state.get(k) for k in ('status','stage','exitCode','elapsedSeconds','error')},ensure_ascii=False))
    return code


if __name__=='__main__':raise SystemExit(main())
