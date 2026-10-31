"""只读完整域已存primitive/守恒量，定位温度/总焓极值，不据单点判定物理资格。"""
from pathlib import Path
import json
import numpy as np
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from su2SSTInterface import state


def main():
    base=ROOT/'outputs/v127_validation';source=Path(read(ROOT/'outputs/v126_validation/main01/engine_folder.json')['folder'])
    report={};locks={}
    for name,folder in [('source49',source),('control149',base/'control'),('cap05_149',base/'cap05')]:
        r=read(folder/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0
        p=folder/'restart_flow.csv';assert sha(p)==r['outputSHA256'][p.name]
        a,cols,nd,rho,u,e=state(p,True,True);c=config(folder/'case.cfg')
        R=float(c['GAS_CONSTANT']);gamma=float(c['GAMMA_VALUE']);Tinf=float(c['FREESTREAM_TEMPERATURE'])
        speed=float(c['MACH_NUMBER'])*np.sqrt(gamma*R*Tinf);h0inf=gamma*R*Tinf/(gamma-1)+speed**2/2+r['physicsModel']['kInfinity']
        pressure=a[:,cols.index('Pressure')];temperature=a[:,cols.index('Temperature')]
        # 使用同阶段p和守恒总能量得到总焓；不错误混入随后更新的k。
        h0=a[:,cols.index('Energy')]/rho+pressure/rho
        fields={'Temperature_K':temperature,'Pressure_Pa':pressure,'TotalEnthalpy_Jkg':h0}
        extrema={}
        for field,v in fields.items():
            ids=np.unique(np.r_[np.argmin(v),np.argsort(v)[-5:]])
            extrema[field]=dict(minimum=float(min(v)),maximum=float(max(v)),samples=[dict(node=int(i),xyzM=a[i,1:4].tolist(),value=float(v[i])) for i in ids])
        report[name]=dict(nodes=len(a),freestreamTotalEnthalpyIncludingK_Jkg=float(h0inf),
            h0Above101ReferenceNodeCount=int(sum(h0>1.01*h0inf)),h0Above110ReferenceNodeCount=int(sum(h0>1.1*h0inf)),
            extrema=extrema,diagnosticOnly=True,physicsQualified=False)
        locks[str(p)]=sha(p);locks[str(folder/'engine_process.json')]=sha(folder/'engine_process.json')
    target=base/'state_extrema.json';assert not target.exists()
    target.write_text(json.dumps(dict(cases=report,sourceSHA256=dict(locks,**{str(Path(__file__).resolve()):sha(Path(__file__).resolve())}),
        interpretation='Nodal diagnostics, not area coverage or an added acceptance gate; no local residual/timestep data inferred.',solverRuns=0),indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:dict(maxT=v['extrema']['Temperature_K']['maximum'],nodesAbove110H0=v['h0Above110ReferenceNodeCount'],
        maxH0=v['extrema']['TotalEnthalpy_Jkg']['maximum'],referenceH0=v['freestreamTotalEnthalpyIncludingK_Jkg']) for k,v in report.items()},indent=2))


if __name__=='__main__':main()
