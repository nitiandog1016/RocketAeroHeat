"""独立逐列核对真实初始化守恒字段与扩散耗散互补条件，不调用初猜求解器。"""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import read,sha
from su2SSTInterface import state


def main():
    source=ROOT/'outputs/v127_validation/control';folder=ROOT/'outputs/v130_validation/column_final'
    target=folder.parent/'initial_state_independent.json';assert not target.exists()
    a,keys,_,rho,_,e=state(source/'restart_flow.csv',True,True)
    b,newkeys,_,_,_,ne=state(folder/'initial_flow.csv',True,True);assert keys==newkeys
    values=np.load(str(folder/'column_values.npz'));ids=values['nodeIDs'];anchors=values['anchorLayers']
    kcol=keys.index('Turb_Kin_Energy');wcol=keys.index('Omega');ecol=keys.index('Energy')
    mask=np.ones(len(a),bool);mask[ids]=False;assert np.array_equal(a[mask],b[mask])
    keep=[j for j in range(len(keys)) if j not in [kcol,wcol,ecol]];assert np.array_equal(a[:,keep],b[:,keep])
    assert np.array_equal(b[ids,ecol],a[ids,ecol]+rho[ids]*(b[ids,kcol]-a[ids,kcol]))
    assert max(abs(ne/e-1))<2e-12 and min(b[ids,kcol])>=1e-10 and min(b[ids,wcol])>1e-4
    werr=0.;kerr=0.;active=0;allids=[];nw=9878
    for j,end in enumerate(anchors):
        if end<2:continue
        ix=np.arange(end+1)*nw+j;allids.extend(ix[1:-1])
        assert np.array_equal(a[ix[[0,-1]]],b[ix[[0,-1]]])
        y=np.linalg.norm(a[ix,1:4]-a[j,1:4],axis=1);k=b[ix,kcol];w=b[ix,wcol]
        nu=a[j,keys.index('Laminar_Viscosity')]/rho[j]
        # 直接在有量纲斜率上复算，与生产代码的无量纲三对角组装独立。
        fluxw=np.diff(w)/np.diff(y);fluxk=np.diff(k)/np.diff(y)
        diffw=nu*2*np.diff(fluxw)/(y[2:]-y[:-2]);diffk=nu*2*np.diff(fluxk)/(y[2:]-y[:-2])
        rw=.075*w[1:-1]**2-diffw;rk=.09*w[1:-1]*k[1:-1]-diffk
        floor=k[1:-1]==1e-10;active+=int(sum(floor))
        wscale=np.maximum(.075*w[1:-1]**2,nu*2*(abs(fluxw[1:])+abs(fluxw[:-1]))/(y[2:]-y[:-2]))
        kscale=np.maximum(.09*w[1:-1]*k[1:-1],nu*2*(abs(fluxk[1:])+abs(fluxk[:-1]))/(y[2:]-y[:-2]))
        werr=max(werr,float(max(abs(rw)/wscale)))
        defect=np.where(floor,np.minimum(rk,0),rk)
        kerr=max(kerr,float(max(abs(defect)/kscale)))
    assert np.array_equal(np.sort(allids),ids)
    assert werr<1e-10 and kerr<1e-10
    assert sha(source/'restart_flow.csv')==sha(folder.parent/'control/initial_flow.csv')
    report=dict(passed=True,selectedNodes=len(ids),unchangedNodes=int(sum(mask)),activeNativeKFloorNodes=active,
        maximumLocalOmegaEquationRelativeError=werr,maximumLocalKComplementarityRelativeDefect=kerr,
        maximumThermalInternalEnergyRelativeDifference=float(max(abs(ne/e-1))),
        controlByteIdentical=True,boundariesAndOtherFieldsUnchanged=True,solverRuns=0,physicsQualified=False,
        sourceSHA256={str(p):sha(p) for p in [source/'restart_flow.csv',folder/'initial_flow.csv',folder/'column_values.npz',Path(__file__).resolve()]})
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
