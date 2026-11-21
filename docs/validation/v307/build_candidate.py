"""在带线热敏珠原第三档体网上只换互不冲突的对角线，四分点测线留给独立审核。"""
import argparse
from pathlib import Path
import shutil

from flip_bound import DEFAULT_BASE,HERE,OUTPUTS,ROOT,START,STOP,LOW,HIGH,read,save,sha,load_mesh,meridian,meridian_faces,wall_frame,np
from thermistorLeadMesh import sweep
from thermistorOuterDomain import write_mesh


def enclosed(xy,tri,wall):
    """距离函数为 1-Lipschitz。覆盖每个原三角形的 16 个子三角形，证明整片都在请求范围内。"""
    def bary(i,j):
        return np.array([4-i-j,i,j])/4
    weights=[]
    for i in range(4):
        for j in range(4-i):
            weights.append((bary(i,j),bary(i+1,j),bary(i,j+1)))
            if i+j<3:weights.append((bary(i+1,j),bary(i+1,j+1),bary(i,j+1)))
    parts=np.einsum("spv,tvd->tspd",np.asarray(weights),xy[tri])
    center=parts.mean(axis=2);radius=np.linalg.norm(parts-center[:,:,None,:],axis=3).max(axis=2).ravel()
    pts=center.reshape(-1,2)
    _,distance,_=wall_frame(pts,wall[START:STOP+1])
    _,aft,_=wall_frame(pts,wall[:START+1])
    _,fore,_=wall_frame(pts,wall[STOP:])
    valid=(distance-radius>LOW)&(distance+radius<HIGH)&(np.minimum(aft,fore)-distance>2*radius)
    return valid.reshape(len(tri),len(weights)).all(axis=1)


def unaltered_diagonal_cuts(points,quads,wall):
    """核对真实边界及两条对角线。射线远离端点且不穿过任一对角线的换边天然符合。"""
    p=points[quads];endpoints=np.array([[0,1],[2,3],[0,2],[2,1],[1,3],[3,0]])
    seg=p[:,endpoints];keep=np.ones(len(p),bool)
    for j in range(START,STOP):
        t=wall[j+1]-wall[j];t/=np.linalg.norm(t)
        n=np.array([-t[1],t[0]]);origin=(wall[j+1]+wall[j])/2
        z=(seg-origin)@t;h=(seg-origin)@n
        den=z[:,:,1]-z[:,:,0];u=np.divide(-z[:,:,0],den,out=np.zeros_like(den),where=abs(den)>1e-14)
        cut=h[:,:,0]+u*(h[:,:,1]-h[:,:,0])
        hit=(abs(den)>1e-14)&(u>=-1e-10)&(u<=1+1e-10)&(cut>LOW+1e-9)&(cut<HIGH-1e-9)
        for k in np.flatnonzero(keep & hit[:,:2].any(axis=1)):
            rim=cut[k,2:][hit[k,2:]]
            def combined(e):
                x=np.sort(np.r_[rim,cut[k,e]] if hit[k,e] else rim)
                return x if len(x)<2 else x[np.r_[True,np.diff(x)>1e-9]]
            old,new=combined(0),combined(1)
            keep[k]=len(old)==len(new) and np.allclose(old,new,atol=1e-11,rtol=0)
    return keep


def select(base,bound):
    summary=read(bound/"bound.json")
    assert summary["meetsUpperGate"] and summary["certificateSHA256"]==sha(bound/"certificate.npz")
    assert sha(base/"mesh.su2")==summary["baseMeshSHA256"] and sha(base/"nodes.npz")==summary["baseMapSHA256"]
    cells,xyz,marks=load_mesh(base/"mesh.su2")
    with np.load(base/"nodes.npz") as z:
        ids=z["ids"];np.testing.assert_array_equal(xyz,z["xyz"])
    xy=xyz[ids[0],:2];tri,_=meridian_faces(marks,ids,len(xyz));wall=meridian()[0]
    with np.load(bound/"certificate.npz") as z:
        np.testing.assert_array_equal(ids[0,tri],z["triNodeIDs"])
        lookup=np.full(len(xyz),-1,int);lookup[ids[0]]=np.arange(len(xy))
        quads=lookup[z["quadNodeIDs"]];pairs=z["pairs"];gain=z["gainFraction"]
        candidates=np.flatnonzero(z["convex"]&z["inside"]&z["rayLockUpperSet"]&(gain>0))
    active=np.unique(pairs[candidates]);valid=np.zeros(len(tri),bool)
    valid[active]=enclosed(xy,tri[active],wall)
    candidates=candidates[valid[pairs[candidates]].all(axis=1)]
    candidates=candidates[unaltered_diagonal_cuts(xy,quads[candidates],wall)]
    used=set();chosen=[];total=0.
    for k in sorted(candidates,key=lambda i:(-gain[i],int(i))):
        owners=set(map(int,pairs[k]))
        if used.isdisjoint(owners):
            chosen.append(k);used.update(owners);total+=gain[k]
            if total>=.20:break
    return dict(cells=cells,xyz=xyz,marks=marks,ids=ids,tri=tri,quads=quads,pairs=pairs,gain=gain,chosen=chosen,
                total=total,strictCandidateCount=len(candidates))


def build(base,bound,out,verify_existing=False):
    model=select(base,bound)
    if not verify_existing:out.mkdir(parents=True,exist_ok=False)
    receipt=out/("selection_verified.json" if verify_existing else "selection.json")
    chosen=np.asarray(model["chosen"],int)
    selection=dict(sourceMeshSHA256=sha(base/"mesh.su2"),certificateSHA256=sha(bound/"certificate.npz"),
                   ruleSHA256=sha(HERE/"CANDIDATE.md"),selectionVerifierSHA256=sha(Path(__file__).resolve()),
                   strictCandidateCount=model["strictCandidateCount"],pairs=chosen.tolist(),predictedFraction=model["total"],
                   constructible=bool(model["total"]>=.20))
    if not selection["constructible"]:
        save(receipt,selection);return
    cells=model["cells"];marks=model["marks"];ids=model["ids"];xyz=model["xyz"];xy=xyz[ids[0],:2]
    cell_id={tuple(sorted(c)):i for i,c in enumerate(cells)}
    face_id={k:{tuple(sorted(f)):i for i,f in enumerate(marks[k])} for k in ("sym0","sym1")}
    assert len(cell_id)==len(cells)
    replaced=[]
    for k in chosen:
        a,b,c,d=model["quads"][k]
        for owner,new in zip(model["pairs"][k],([c,d,a],[d,c,b])):
            old=model["tri"][owner];points=xy[new]
            if np.cross(points[1]-points[0],points[2]-points[0])<0:new[1],new[2]=new[2],new[1]
            new=np.asarray(new,int)
            for low,high in zip(ids,ids[1:]):
                index=cell_id[tuple(sorted(np.r_[low[old],high[old]]))]
                cells[index]=tuple(sweep(new,low,high));replaced.append(index)
            for tag,ring,oriented in (("sym0",ids[0],new[[0,2,1]]),("sym1",ids[-1],new)):
                index=face_id[tag][tuple(sorted(ring[old]))]
                marks[tag][index]=tuple(ring[oriented])
    assert len(set(replaced))==len(replaced)==8*len(chosen)
    if verify_existing:
        actual,coords,markers=load_mesh(out/"mesh.su2")
        assert cells==actual and marks==markers and sha(out/"nodes.npz")==sha(base/"nodes.npz")
        np.testing.assert_array_equal(coords,xyz)
    else:
        write_mesh(out/"mesh.su2",cells,xyz,marks)
        shutil.copyfile(base/"nodes.npz",out/"nodes.npz")
    selection.update(replacedCells=len(replaced),candidateMeshSHA256=sha(out/"mesh.su2"))
    save(receipt,selection);print(len(chosen),model["total"],len(replaced))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("bound",type=Path);p.add_argument("output",type=Path)
    p.add_argument("--base",type=Path,default=DEFAULT_BASE)
    p.add_argument("--verify-existing",action="store_true",help="只校核已写出的全部连接并另存回执");a=p.parse_args()
    bound=a.bound.resolve();out=a.output.resolve();base=a.base.resolve()
    for f in (out,bound):f.relative_to(OUTPUTS)
    base.relative_to(ROOT/"outputs");build(base,bound,out,a.verify_existing)
