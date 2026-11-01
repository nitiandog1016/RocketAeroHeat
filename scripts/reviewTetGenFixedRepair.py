"""只读独立复核TetGen小例失败：原生数组、固定面、体积与逆矩阵质量。"""
from pathlib import Path
import json,time
import numpy as np
from prepareTetGenRuntime import ROOT,DEST,sha
from optimizeSU2FanExterior import sicn

def read_raw(path):
    lines=path.read_text().splitlines();n,t,f=map(int,lines[0].split())
    assert len(lines)==1+n+t+f
    points=np.array([list(map(float,s.split())) for s in lines[1:1+n]])
    cells=np.array([list(map(int,s.split())) for s in lines[1+n:1+n+t]],np.int64)
    faces=np.array([list(map(int,s.split())) for s in lines[1+n+t:]],np.int64)
    return points,cells,faces

def keys(a):return sorted(tuple(sorted(row)) for row in a.tolist())

def independent_topology(P,T):
    groups={}
    for t in T:
        for idx in [(0,2,1),(0,1,3),(1,2,3),(2,0,3)]:
            f=tuple(int(t[j]) for j in idx);key=tuple(sorted(f))
            k=f.index(min(f));direction=1 if f[(k+1)%3]>f[(k-1)%3] else -1
            groups.setdefault(key,[]).append(direction)
    boundary=[]
    for face,sign in groups.items():
        assert len(sign) in [1,2]
        if len(sign)==2:assert sum(sign)==0
        else:boundary.append(face)
    volume=float(np.linalg.det(np.stack([P[T[:,1]]-P[T[:,0]],P[T[:,2]]-P[T[:,0]],P[T[:,3]]-P[T[:,0]]],axis=2)).sum()/6)
    return boundary,volume

def main():
    start=time.monotonic();folder=ROOT/'docs/validation/v142/hollow_test';report=json.loads((folder/'report.json').read_text())
    runtime=json.loads((DEST/'runtime.json').read_text());native=folder/'native';receipt=json.loads((native/'process/process.json').read_text())
    paths=[Path(__file__).resolve(),ROOT/'scripts/optimizeSU2FanExterior.py',folder/'report.json',folder/'fixture.npz',
        native/'input.dat',native/'native_import.dat',native/'native_result.dat',native/'native_options.json',
        native/'process/process.json',native/'process/stdout.log',DEST/'runtime.json']
    locked={str(p.resolve()):sha(p) for p in paths}
    assert receipt['exitCode']==0 and all(receipt['sourceStableAfter'].values())
    for r in [receipt,report]+runtime['receipts']:
        for p,h in r['sourceSHA256'].items():assert sha(Path(p))==h
    assert sha(Path(runtime['executable']))==runtime['executableSHA256']
    for n,r in json.loads((DEST/'source.json').read_text())['sourceGitBlobChecks'].items():assert sha(DEST/'source'/n)==r['sha256']
    P,T,F=read_raw(native/'input.dat');I,J,K=read_raw(native/'native_import.dat');Q,U,G=read_raw(native/'native_result.dat')
    assert np.array_equal(P,I) and np.array_equal(T,J) and np.array_equal(F,K)
    fixture=np.load(str(folder/'fixture.npz'));assert np.array_equal(P[:,:3],fixture['vertices']) and np.array_equal(T,fixture['tetrahedra'])
    assert np.array_equal(P,Q) and keys(T)==keys(U) and keys(F[:,:3])==keys(G[:,:3])
    assert np.all(F[:,3]==1) and np.all(G[:,3]==1)
    before,vol0=independent_topology(P[:,:3],T);after,vol1=independent_topology(Q[:,:3],U)
    assert sorted(before)==keys(F[:,:3]) and sorted(after)==keys(G[:,:3])
    assert abs(vol0-26)<1e-12 and abs(vol1-26)<1e-12
    boundary=np.unique(F[:,:3]);assert np.array_equal(P[boundary,3],boundary+1)
    assert np.array_equal(np.flatnonzero(P[:,3]==0),[64])
    q0=sicn(P[:,:3],T);q1=sicn(Q[:,:3],U);bad=np.flatnonzero(q1<.1)
    err=abs(float(min(q1))-report['minimumSICNAfter']);assert err<1e-12
    assert not report['passed'] and len(bad)==1 and min(q1)==min(q0)
    nativeOptions=json.loads((native/'native_options.json').read_text());assert nativeOptions['switches']=='rY1q2/5O2/7S10000zVC'
    assert nativeOptions['optLevel']==2 and nativeOptions['optScheme']==7 and nativeOptions['maxDihedral']==165
    text=(native/'process/stdout.log').read_text();assert 'Added 0 (0,0,0) Steiner points.' in text
    assert not (ROOT/'outputs/v142_validation').exists(),'unexpected real candidate outputs'
    result=dict(passed=True,reviewScope='independent verification of rejection, not meshing qualification',
        exactRawImport=True,allNativePointsAndMarkersExactlyUnchanged=True,allNativeTetConnectivityUnchanged=True,
        fixedBoundaryFacesExactlyUnchanged=True,originalInteriorPointCount=1,additionalPoints=0,
        beforeVolume=vol0,afterVolume=vol1,knownCavityVolume=1,originalBoundaryFaces=len(F),
        independentMinimumSICN=float(min(q1)),independentSICNMaximumMinimumDifference=err,
        badTetCount=len(bad),badNativeTetRowsOneBased=(bad+1).tolist(),badNativeTetNodeIDsOneBased=(U[bad]+1).tolist(),
        nativeExitCode=0,testExitCode=1,realCandidateExecuted=False,nativeRunCount=1,
        noNewMATLABOrSU2OrUITests=True,sourceSHA256=locked,
        sourceStableAfter={p:sha(Path(p))==h for p,h in locked.items()},seconds=time.monotonic()-start,
        physicsQualified=False,newQualifiedArea=0)
    assert all(result['sourceStableAfter'].values())
    dst=folder/'independent_review.json';assert not dst.exists();dst.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['sourceSHA256','sourceStableAfter']},indent=2))

if __name__=='__main__':main()
