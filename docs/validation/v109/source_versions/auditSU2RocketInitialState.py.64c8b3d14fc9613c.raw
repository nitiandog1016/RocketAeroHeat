"""从CSV和几何级数独立复核三维初始猜测，不能将初值当边界层解。"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder');folder=Path(parser.parse_args().folder).resolve()
    report=json.loads((folder/'initialization_manifest.json').read_text())
    for p,h in report['inputSHA256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    cfg=json.loads((folder/'case_manifest.json').read_text());f=cfg['flight']
    state=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1)
    rho=state[:,4];velocity=state[:,5:8]/rho[:,None];internal=state[:,8]/rho-(velocity**2).sum(axis=1)/2
    T=internal*.4/287.05287;pressure=rho*287.05287*T
    # 独立采用等比层距公式，不重复生成器按三维位移求比例的方法。
    eta=(1.25**np.arange(25)-1)/(1.25**24-1);expected=np.r_[eta**2*(3-2*eta),np.ones(16)]
    expected=np.repeat(expected,report['wallNodes']);u=expected[:,None]*cfg['expectedVelocity']
    t=f['wallTemperature']+(f['temperature']-f['wallTemperature'])*expected
    errors=dict(velocity=float(np.max(abs(velocity-u))/np.linalg.norm(cfg['expectedVelocity'])),
                temperature=float(np.max(abs(T-t))/f['wallTemperature']),pressure=float(np.max(abs(pressure/f['pressure']-1))))
    assert max(errors.values())<2e-10 and np.min(rho)>0 and np.min(internal)>0
    assert np.all(velocity[:report['wallNodes']]==0)
    a=dict(independentReference='geometric layer ratio and conservative-state EOS',errors=errors,
           positiveInternalEnergy=True,noSlip=True,physicsQualified=False,
           sourceSHA256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                         (Path(__file__).resolve(),folder/'initial_flow.csv',folder/'initialization_manifest.json')})
    (folder/'initial_audit.json').write_bytes((json.dumps(a,indent=2)+'\n').encode('utf8'));print(json.dumps(errors))


if __name__=='__main__':main()
