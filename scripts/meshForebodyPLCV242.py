"""模块职责：唯一成熟TetGen调用，先保存全部原生结果再作几何与原质量门槛审核。"""
import itertools,subprocess
import numpy as np
from scipy.sparse.csgraph import connected_components
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from tetgenFixedInterface import read_native,tet_audit,sorted_cells
from cgalTetrahedralContract import faces_and_neighbors
from reviewSU2ExplicitQuality import quality
from clipEulerForebody import write_su2
BASE=ROOT/'outputs/v242_validation'


def mesh():
    out=BASE/'tetgen01';out.mkdir(exist_ok=False);runtime=read(BASE/'build01/runtime.json');assert runtime['passed']
    exe=runtime['executable'];assert sha(exe)==runtime['executableSHA256']
    z=np.load(str(BASE/'surface01/plc.npz'));p=z['vertices'];f=z['faces'];markers=z['markers'];inp=BASE/'surface01/plc.dat'
    assert sha(inp)==read(BASE/'surface01/report.json')['plcSHA256']
    cmd=[exe,str(inp),str(out/'native')];save(out/'command.json',cmd)
    result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(out/'tetgen.txt').write_bytes(result.stdout)
    assert result.returncode==0,('TETGEN_FAILED',result.returncode)
    # 实际导入数组必须同输入逐位相同；不以输入文件自身回读代替真实核心入口。
    with (out/'native_import.dat').open() as stream:
        nv,nf,nh=map(int,stream.readline().split());assert (nv,nf,nh)==(len(p),len(f),0)
        ip=np.loadtxt(itertools.islice(stream,nv));iff=np.loadtxt(itertools.islice(stream,nf),dtype=np.int64);assert not stream.read().strip()
    assert np.array_equal(ip[:,:3],p) and np.array_equal(ip[:,3],np.arange(1,len(p)+1))
    assert np.array_equal(iff[:,:3],f) and np.array_equal(iff[:,3],markers)
    q,t,g,ref,mark=read_native(out/'native_result.dat');boundary,geo=tet_audit(q,t)
    assert np.array_equal(sorted_cells(boundary),sorted_cells(g));fixed=ref>0
    assert np.array_equal(np.sort(ref[fixed]),np.arange(1,len(p)+1)) and np.all(fixed[np.unique(g)])
    error=float(np.max(abs(q[fixed]-p[ref[fixed]-1])));tol=16*np.finfo(float).eps*max(1.,float(np.ptp(p,axis=0).max()));assert error<=tol
    node_map=np.empty(len(q),np.int64);node_map[fixed]=ref[fixed]-1;node_map[~fixed]=np.arange(len(p),len(p)+sum(~fixed))
    r=np.r_[p,q[~fixed]];u=node_map[t];newf=node_map[g]
    assert np.array_equal(sorted_cells(newf),sorted_cells(f))
    expected={tuple(sorted(face)):int(m) for face,m in zip(f,markers)}
    assert all(expected[tuple(sorted(face))]==m for face,m in zip(newf,mark))
    _,graph=faces_and_neighbors(u);assert connected_components(graph,directed=False,return_labels=False)==1
    volume_error=abs(geo['volume']/float(z['expectedVolume'])-1);assert volume_error<2e-10
    qvalues=quality(r[u]);good=bool(qvalues.min()>=.1);sizeok=len(u)<=200000 and len(r)<=50000;assert sizeok
    np.savez_compressed(str(out/'mesh.npz'),vertices=r,tetrahedra=u,wall=z['wall'],outlet=f[markers==2],farfield=f[markers==3],
      originalFaceIDs=z['originalFaceIDs'],originalArea=z['originalArea'],expectedVolume=z['expectedVolume'],SICN=qvalues)
    write_su2(out/'mesh.su2',dict(xyz=r,cells=u,markers=dict(wall=z['wall'],outlet=f[markers==2],farfield=f[markers==3])))
    save(out/'report.json',dict(passed=True,geometryPassed=True,qualityPassed=good,sizePassed=sizeok,nodes=len(r),cells=len(u),
      nativeZeroHoleInputExact=True,nativeBoundaryMarkersExact=True,originalWallExactlyPreserved=True,
      minimumSICN=float(qvalues.min()),cellsBelow01=int(sum(qvalues<.1)),qualityPercentiles=np.percentile(qvalues,[0,1,50,99,100]).tolist(),
      maximumNativeFixedNodeDifferenceM=error,roundtripToleranceM=tol,originalCoordinatesRestoredOnlyWithinTolerance=True,
      volumeRelativeDifference=volume_error,nativeGeometry=geo,netgenNeeded=not good,physicalQualification=False,newQualifiedAreaM2=0,
      inputSHA256=sha(inp),runtimeSHA256=sha(BASE/'build01/runtime.json'),meshSHA256=sha(out/'mesh.su2')))


if __name__=='__main__':raise SystemExit('v242 fixed task closed; use NEXT_DECISION')
