# -*- coding: utf-8 -*-
"""读取冻结C与真实质心，独立审计给定流向排序及上下三角五量块权重。"""
from pathlib import Path
import os,json,hashlib
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'flow_coupling_diagnostic01.json';assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
geometry=VAL/'flowordergeometry01/flow_ordering_geometry.mat';g=loadmat(geometry,simplify_cells=True);center=g['centroid'];direction=g['unitFlowDirection'];volume=g['volume'];report=g['report']
source=Path(g['sourceFile']);assert sha(source)==report['sourceSHA256'];s=loadmat(source,simplify_cells=True,variable_names=['domain','fs'])
np.testing.assert_array_equal(center,s['domain']['cells']['centroid']);np.testing.assert_array_equal(volume,s['domain']['cells']['volume']);np.testing.assert_array_equal(direction,s['fs']['unitFlowDirection']);del s
projection=center@(direction/np.linalg.norm(direction));n=len(center);flow_order=np.lexsort((np.arange(n),projection));inverse=np.empty(n,int);inverse[flow_order]=np.arange(n)
packedfile=VAL/'blockilureal01/steady_graph_packed_evidence.mat';p=loadmat(packedfile,simplify_cells=True,variable_names=['packed'])['packed'];order=p['order'].astype(int)-1;ptr=p['rowPointer'].astype(np.int64)-1;col=p['column'].astype(np.int64)-1;owners=np.repeat(np.arange(n),np.diff(ptr))
assert n==p['cellCount']==report['cellCount'] and p['sourceFaces']==report['sourceFaces'] and p['layers']==report['layers']
metrics={name:{k:np.zeros(5) for k in ('diagonal','lower','upper','offDiagonalProjectionTie')} for name in ('symrcm','givenFlowProjection')};counts={name:{k:0 for k in metrics[name]} for name in metrics}
for first in range(0,len(col),4096):
 ids=slice(first,first+4096);rr=owners[ids];cc=col[ids];physical_r=order[rr];physical_c=order[cc];weights=np.sum(p['values'][:,:,ids]**2,axis=1).T
 flow_r=inverse[physical_r];flow_c=inverse[physical_c];ties=(physical_r!=physical_c)&(projection[physical_r]==projection[physical_c])
 for name,r,c in [('symrcm',rr,cc),('givenFlowProjection',flow_r,flow_c)]:
  for key,mask in [('diagonal',r==c),('lower',r>c),('upper',r<c),('offDiagonalProjectionTie',ties)]:metrics[name][key]+=weights[mask].sum(axis=0);counts[name][key]+=int(np.count_nonzero(mask))
outmetrics={}
for name,m in metrics.items():
 denominator=m['diagonal']+m['lower']+m['upper'];off=m['lower']+m['upper']
 outmetrics[name]=dict(blockCounts=counts[name],blockFrobeniusSquaredFraction={k:float(v.sum()/denominator.sum()) for k,v in m.items()},
  lowerFractionOfOffDiagonalByConservationRow=(m['lower']/np.maximum(off,np.finfo(float).tiny)).tolist(),
  lowerFractionOfAllOffDiagonal=float(m['lower'].sum()/off.sum()))
out=dict(status='passed-source-domain-centroid-identity-and-complete-C-projected-order-coupling-diagnostic',processID=os.getpid(),cellCount=n,
 givenUnitFlowDirection=direction.tolist(),projectionRangeM=[float(projection.min()),float(projection.max())],metrics=outmetrics,
 interpretation='These are squared normalized-C block weights, not eigenvalues or pressure/heat errors. Flow projection uses given freestream, not local streamlines. Ordering usefulness requires actual factor/action tests.',
 newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
 inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (source,geometry,packedfile)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
