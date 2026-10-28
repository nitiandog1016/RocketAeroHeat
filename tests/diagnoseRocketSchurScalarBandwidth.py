# -*- coding: utf-8 -*-
"""用原图超图与非零同CV系数复算等价标量排列的初始带宽界。"""
from pathlib import Path
import os,json,hashlib
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'schur_scalar_bandwidth02.json';assert not target.exists()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
file=VAL/'blockfillreal01/expanded_pattern_evidence.mat';s=loadmat(file,simplify_cells=True,variable_names=['pattern']);p=s['pattern'];n=len(p['order']);owners=np.repeat(np.arange(n),np.diff(p['rowPointer'].astype(np.int64)));col=p['column'].astype(np.int64)-1;b=int(np.max(abs(col-owners)))
proofFile=VAL/'schurcheckpointreal01_independent_real_checkpoint.json';proof=json.loads(proofFile.read_text());assert proof['originalSchurDiagonal'][0][1]!=0 and p['enrichmentReport']['allOriginalBlockValuesKept']
out=dict(status='source-bound-reproduced-structural-bandwidth-bounds-for-equivalent-Schur-scalar-orderings',processID=os.getpid(),cellCount=n,
 originalPatternIsSubsetOfAuditedLevel1=True,level1BlockBandwidth=b,fullTwoHopSchurBlockBandwidthUpperBound=2*b,
 fieldMajorScalarBandwidthLowerBound=n,cellInterleavedScalarBandwidthUpperBound=4*b+1,lowerBoundOnInitialBandwidthReduction=n/(4*b+1),
 numericalAccuracyFillOrRuntimeImprovementClaimed=False,newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,
 inputSHA256={str(f.relative_to(ROOT)):sha(f) for f in (file,proofFile)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
