"""模块职责：只读核对真实重启、全局节点映射与能量合同拒绝，不重跑控制。"""
import csv,copy,json
import numpy as np
from runSU2CompatibilityV221 import ROOT,BASE,BUILD,read,save,sha
from su2PhaseAirContract import read_restart,bind
def main():
    out=BASE/'restart_audit01';out.mkdir(exist_ok=False);origin=BASE/'control01';restart=BASE/'restart01'
    assert read(BASE/'process_restart/report.json')['passed']
    def load(p):
        with p.open(encoding='utf8') as f:return list(csv.DictReader(f))
    first=load(origin/'saved.csv');second=load(restart/'restarted.csv')
    def mapped(rows):
        result={tuple(float(r[k]) for k in ('x','y','z')):np.array([float(r[k]) for k in ('rho','rhou','rhov','rhow','rhoE','nutilde')]) for r in rows}
        assert len(result)==len(rows)==24;return result
    a,b=mapped(first),mapped(second);assert set(a)==set(b)
    for xyz in a:assert a[xyz].tobytes()==b[xyz].tobytes(),xyz
    cfg=(origin/'case.cfg').read_text(encoding='utf8');assert (restart/'case.cfg').read_text(encoding='utf8')==cfg.replace('RESTART_SOL= NO','RESTART_SOL= YES')
    metadata=read(origin/'restart_flow.rah.json')
    expected=bind(origin/'restart_flow.csv',origin/'mesh.su2',origin/'case.cfg',ROOT/'outputs/v220_validation/build02/SU2_CFD.exe',BUILD/'compatibility.exe','manufactured-control')
    frame=read_restart(restart/'restart_flow.csv',metadata,expected,restart/'mesh.su2')
    assert frame['wallHeatFlux'] is None and frame['solidTemperature'] is None
    # 对真实文件做独立小副本：即使更新了载荷SHA，也不能漏点/重复编号/错坐标。
    raw=(origin/'restart_flow.csv').read_text(encoding='utf8').splitlines(True);rejects=[]
    altered={
      'missing_node':raw[:-1],
      'duplicate_node':raw[:-1]+[raw[1]],
      'wrong_coordinate':raw[:1]+[raw[1].replace('1.0000000000000000e+00','1.0001000000000000e+00',1)]+raw[2:]}
    for name,lines in altered.items():
        p=out/(name+'.csv');p.write_text(''.join(lines),encoding='utf8')
        expected2=copy.deepcopy(expected);expected2['payloadSHA256']=sha(p)
        try:read_restart(p,expected2,expected2,restart/'mesh.su2')
        except ValueError as exc:rejects.append(dict(case=name,reason=str(exc)))
        else:raise AssertionError(name)
    save(out/'report.json',dict(passed=True,actualRestartNodesBitwise=24,conservativePlusSAValuesBitwise=144,
        globalIDCoordinateCoverage=True,localOrderIsNotGlobalOrder=True,uncomputedDiagnosticsDiscarded=True,
        rejectionCases=rejects,boundaryAndSAScope=read(BASE/'audit01/report.json'),
        interfaceControlAccepted=True,manufacturedControlOnly=True,CFDIterations=0,newQualifiedAreaM2=0,
        fullFlowJacobianVerified=False,interiorSADiffusionIndependentlyVerified=False))
    print(json.dumps(dict(passed=True,actualRestartNodesBitwise=24,additionalRejections=len(rejects))))
if __name__=='__main__':main()
