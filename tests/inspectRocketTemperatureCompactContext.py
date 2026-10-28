"""读取已保存的当前QR/全部H，独立检查完整组支撑；不调用NS或借用旧状态。"""
from pathlib import Path
import os, sys, json
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import sha

file = Path(sys.argv[1]).resolve(); target = file.parent/'independent_compact_context.json'; assert not target.exists()
before = sha(file); saved = loadmat(file, simplify_cells=True, variable_names=['factor','rows','preparationSHA256'])
factor=saved['factor']; rows=saved['rows']; H=rows['H'].tocsr(); n=int(factor['cellCount'])
representatives=[]; mapping=[]
for i in range(H.shape[0]):
    row=H.getrow(i); match=next((k for k,r in enumerate(representatives) if (row-H.getrow(r)).nnz==0),None)
    if match is None: match=len(representatives); representatives.append(i)
    mapping.append(match)
touched=np.unique(H.indices % n); covered=np.zeros(n,dtype=bool); plans=[]
buckets=factor['buckets']; buckets=[buckets] if isinstance(buckets,dict) else list(np.atleast_1d(buckets))
for k,b in enumerate(buckets,1):
    count=np.size(b['patchIDs'])
    if not count: continue
    members=np.asarray(b['members'],int).reshape(k,count,order='F')-1
    pages=np.flatnonzero(np.any(np.isin(members,touched),axis=0))
    maps=np.asarray(b['map']).reshape(5*k,5*k,count,order='F')
    for page in pages:
        cv=members[:,page]; assert not np.any(covered[cv]); covered[cv]=True
        ids=(cv[:,None]+n*np.arange(5)).ravel(); plans.append((ids,maps[:,:,page]))
assert np.all(covered[touched]); cv=np.flatnonzero(covered)
ids=(cv[:,None]+n*np.arange(5)).ravel(order='F'); lookup={value:i for i,value in enumerate(ids)}
columns=np.zeros((len(ids),len(representatives)))
for original,P in plans:
    compact=np.array([lookup[x] for x in original]); columns[compact]=P.T@H[representatives][:,original].toarray().T
norms=np.linalg.norm(columns,axis=0); assert np.all(norms>0)
U,R=np.linalg.qr(columns/norms,mode='reduced'); pivots=abs(np.diag(R)); assert np.all(pivots>1e-10)
assert sha(file)==before
out=dict(status='passed-independent-current-complete-QR-group-support-and-all-original-row-span',
    inputFile=str(file),inputSHA256=before,preparationSHA256=saved['preparationSHA256'],
    originalRows=H.shape[0],exactDistinctRows=len(representatives),originalToUniqueRow=(np.array(mapping)+1).tolist(),
    sourceSupportCVs=len(touched),completeQRSupportCVs=len(cv),storedCoordinates=len(ids),
    basisArrayBytes=U.nbytes,equivalentFullBasisArrayBytes=8*5*n*len(representatives),
    normalizedQRDiagonal=pivots.tolist(),orthogonalityError=float(np.linalg.norm(U.T@U-np.eye(len(representatives)),ord=np.inf)),
    newActualNSEvaluations=0,newFlowAcceptances=0,physicalSolution=False,newProductionSolvedAreaM2=0,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8'); print(json.dumps(out))
