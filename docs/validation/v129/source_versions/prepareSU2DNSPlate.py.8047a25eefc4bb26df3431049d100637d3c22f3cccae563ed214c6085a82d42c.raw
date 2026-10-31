"""从MATLAB生成的平板准备有明确入口假设的M6冷壁SST算例。"""
from pathlib import Path
import argparse,json
from prepareSU2SSTCase import prepare
from compareSU2DNSPlate import reference,REFERENCE_SHA,ROOT
from prepareSU2CheckpointContinuation import config,read,sha


def prepare_case(folder):
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    ref,_=reference(); c=config(folder/'case.cfg')
    for key,val in [('GAS_CONSTANT',287),('PRANDTL_LAM',.71),('FREESTREAM_TEMPERATURE',55.2)]:
        assert float(c[key])==val
    assert abs(float(c['MACH_NUMBER'])/ref['MachFromSpeed']-1)<1e-14
    assert abs(float(c['FREESTREAM_PRESSURE'])/ref['pressurePa']-1)<1e-14
    prepare(folder,intensity=.01,viscosity_ratio=10,prt=.9,
            provenance='M6Tw025 matched-ReTheta developing plate; uniform I=1%, mu_t/mu=10 are research assumptions, not DNS inlet measurements; fully turbulent, no transition')
    p=folder/'case_manifest.json';m=read(p)
    m['publicBenchmark']=dict(referenceSHA256=REFERENCE_SHA,reference=ref,
        matchingRule='unique-ReTheta-2121-interior-v1',exactDNSReplication=False,
        comparisonClass='matched-ReTheta-developing-flat-plate-with-inlet-uncertainty',
        initialState='native uniform freestream; no downloaded flow initialization',
        iterationBudget=10000,nativeSecondsBudget=720,threads=2,
        planSHA256=sha(ROOT/'docs/validation/v129/PLAN.md'))
    p.write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return m['publicBenchmark']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',required=True);a=p.parse_args()
    print(json.dumps(prepare_case(a.folder),ensure_ascii=False))
