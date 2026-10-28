# -*- coding: utf-8 -*-
"""独立核对MAT源到真实原面JSON，并绘制科研PNG；不是浏览器或WebGL视觉验收。"""
from pathlib import Path
import hashlib, json, os, sys
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
root=Path(__file__).resolve().parents[1];val=root/'outputs/v82_validation'
folder=val/(sys.argv[2] if len(sys.argv)>2 else 'workbench01')
suffix=sys.argv[1] if len(sys.argv)>1 else '';assert not suffix or suffix.isdigit()
target=folder/('independent_original_face_audit'+suffix+'.json');assert not target.exists()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
domain_file=root/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
bindings={str(domain_file.relative_to(root)):sha(domain_file)}
d=loadmat(domain_file,simplify_cells=True,variable_names=['domain'])['domain'];mesh=d['inputGeometry'];w=d['walls']
assert np.array_equal(w['originalFaceID'],np.arange(1,19753))
payload_file=folder/'exact-study-data.json';bindings[str(payload_file.relative_to(root))]=sha(payload_file)
payload=json.loads(payload_file.read_text(encoding='utf8'));data=[];checks=[]
for index,c in enumerate(payload['cases']):
    file=folder/c['resultsJSONURL'];bindings[str(file.relative_to(root))]=sha(file)
    view=json.loads(file.read_text(encoding='utf8'));data.append(view)
    np.testing.assert_array_equal(view['vertices'],mesh['V']);np.testing.assert_array_equal(view['faces'],mesh['F']-1)
    if index==0: source=root/'outputs/v81_validation/arithmetic02/steady_krylov_evidence.mat'
    elif index==1: source=val/'integrate01/exact_integration_evidence.mat'
    else: source=val/c['summary'].get('resultEvidenceFolder','exact_krylov01')/'steady_krylov_evidence.mat'
    bindings[str(source.relative_to(root))]=sha(source)
    s=loadmat(source,simplify_cells=True);r=s.get('result',s);info=r['info'];rate=r['rate'];gas=s['gas']
    P=info['primitive'];rho=P[:,0];u=P[:,1:4];cv=gas['R']/(gas['gamma']-1);owner=w['ownerCell'].astype(int)-1
    temperature_rate=(rate[:,4]-np.sum(u*rate[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*P[:,4])*rate[:,0])/(rho*cv)
    expected=dict(pressure=info['wallPressure']/1000,heatFlux=info['wallHeatFluxIntoWall']/1000,
                  wallAdjacentTemperature=P[owner,4],temperatureRate=temperature_rate[owner],wallInput=info['wallTemperatureInput'])
    fields={x['key']:x for x in view['fields']};errors={}
    for key,values in expected.items():
        actual=np.asarray(fields[key]['values']);error=float(np.max(abs(actual-values))/max(1,np.max(abs(values))))
        assert error<3e-13;errors[key]=error
    for key in ['edgeTemperature','recoveryTemperature','wallSolved']:
        assert all(x is None for x in fields[key]['values']) and set(fields[key]['status'])=={0}
    assert not view['researchIdentity']['productionQualified'] and not view['researchIdentity']['HLLCHistoryAppended']
    checks.append(dict(case=c['id'],sourceRelativeDifferences=errors))
    del s,r,info,rate

# 模块二：同一原面与共同全范围色标，保持真实几何比例；PNG只作独立科研图证据。
V=np.asarray(data[0]['vertices']);F=np.asarray(data[0]['faces'],dtype=int)
view=np.asarray([.15,-1.,.4]);view/=np.linalg.norm(view)
right=np.cross([0.,0.,1.],view);right/=np.linalg.norm(right);up=np.cross(view,right)
projected=np.column_stack((V@right,V@up));depth=V@view;order=np.argsort(depth[F].mean(axis=1))
triangles=projected[F][order]
fig=plt.figure(figsize=(14,5.4),layout='constrained')
for row,(key,unit) in enumerate([('pressure','kPa'),('heatFlux','kW/m2')]):
    normalizer=Normalize(*payload['commonRanges'][key]);axes=[]
    for col,index in enumerate([0,2]):
        ax=fig.add_subplot(2,2,row*2+col+1);axes.append(ax)
        fields={x['key']:x for x in data[index]['fields']};values=np.asarray(fields[key]['values'])
        poly=PolyCollection(triangles,facecolors=plt.colormaps['viridis'](normalizer(values[order])),edgecolors='none',rasterized=True)
        ax.add_collection(poly);low=projected.min(axis=0);high=projected.max(axis=0);pad=.03*(high-low)
        ax.set(xlim=(low[0]-pad[0],high[0]+pad[0]),ylim=(low[1]-pad[1],high[1]+pad[1]),xlabel='Orthographic horizontal coordinate (m)',ylabel='View height (m)')
        ax.set_aspect('equal',adjustable='box');ax.tick_params(labelsize=8)
        path_label='primitive tangent' if payload['cases'][2]['summary']['lineage'].get('nonlinearTrialStateCoordinates')=='primitive-tangent' else 'conservative'
        ax.set_title(('Frozen HLLC source' if index==0 else 'Exact Godunov: '+path_label+', 3 corrections')+' / '+key,fontsize=10)
    fig.colorbar(ScalarMappable(norm=normalizer,cmap='viridis'),ax=axes,shrink=.65,label=key+' ('+unit+'); common full range')
fig.suptitle('M5 / 30 km / alpha=20 deg / prescribed Tw=300 K / UNCONVERGED RESEARCH\nOrthographic 3D surface; original geometry proportions preserved',fontsize=12)
comparison_name='original_face_comparison'+suffix+'.png'
residual_name='new_equation_residuals'+suffix+'.png'
fig.savefig(folder/comparison_name,dpi=150);plt.close(fig)
history=payload['cases'][2]['history'];fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
ax.semilogy(history['iteration'],history['maximumResidualPerSecond'],'o-',label='Maximum scaled residual')
ax.semilogy(history['iteration'],history['rmsMerit'],'s-',label='RMS merit')
ax.set(xlabel='Accepted corrections of the NEW equation (not physical time)',ylabel='Residual (1/s)',title='Exact Godunov + full 3D NS / not steady')
ax.grid(True,alpha=.25);ax.legend();fig.savefig(folder/residual_name,dpi=150);plt.close(fig)
assert all(sha(root/p)==h for p,h in bindings.items())
report=dict(status='passed-independent-MAT-to-original-facet-all-field-identity',processID=os.getpid(),cases=checks,
            sourceFaces=19752,scientificPNGs=[comparison_name,residual_name],
            browserVisualOrWebGLRuntimeVerified=False,physicalSolution=False,newProductionSolvedAreaM2=0,
            inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8');print(json.dumps(report))
