"""独立积分棱柱多面体与面配对，检查完整域；不把拓扑合格当作物理验证。"""
import numpy as np

def polyvolume(P,F):
    return float(np.sum(np.einsum('ij,ij->i',P[F[:,0]],np.cross(P[F[:,1]],P[F[:,2]])))/6)

def paired(faces,boundary):
    """从最小节点的前后邻点核对循环方向，独立于单元生成顺序。"""
    key,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
    k=np.argmin(faces,axis=1);row=np.arange(len(faces));n=faces.shape[1]
    signs=np.sign(faces[row,(k+1)%n]-faces[row,(k-1)%n])
    direction=np.bincount(inv,weights=signs);inner=count==2
    assert np.all(count<=2) and np.all(direction[inner]==0),'非流形或内面方向不匹配'
    b=np.sort(boundary,axis=1);order=np.lexsort(b.T[::-1])
    assert np.array_equal(key[count==1],b[order]),'声明边界与体单元边界不同'
    return dict(internalFaces=int(sum(inner)),boundaryFaces=int(sum(count==1)),oppositeOrientation=True)

def audit(data):
    P=data['vertices'];W=data['wedges'].astype(int)-1;T=data['tetrahedra'].astype(int)-1
    F=data['sourceFaces'].astype(int)-1;far=data['farfield'].astype(int)-1
    V=data['sourceVertices'];assert np.array_equal(P[:len(V)],V)
    ref=P[W].mean(axis=1);volume=np.zeros(len(W));closure=np.zeros_like(ref);scale=np.zeros(len(W))
    for face in [[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]]:
        ids=W[:,face];center=P[ids].mean(axis=1)-ref
        for j in range(len(face)):
            a=P[ids[:,j]]-ref;b=P[ids[:,(j+1)%len(face)]]-ref
            av=np.cross(b-a,center-a)/2
            volume+=np.einsum('ij,ij->i',a,av)/3;closure+=av;scale+=np.linalg.norm(av,axis=1)
    assert np.all(volume>0)
    closureError=float(max(np.linalg.norm(closure,axis=1)/scale));assert closureError<2e-12
    q=P[T];tv=np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0]))/6
    assert np.all(tv>0)
    tri=np.concatenate([W[:,[0,2,1]],W[:,[3,4,5]],T[:,[0,2,1]],T[:,[0,1,3]],T[:,[1,2,3]],T[:,[2,0,3]]])
    quad=np.concatenate([W[:,[0,1,4,3]],W[:,[1,2,5,4]],W[:,[2,0,3,5]]])
    tr=paired(tri,np.concatenate([F[:,[0,2,1]],far]));qr=paired(quad,np.empty((0,4),int))
    outer=polyvolume(P,far);body=polyvolume(V,F);total=float(sum(volume)+sum(tv))
    relative=abs(total-(outer-body))/total;assert relative<2e-10
    return dict(passed=True,prisms=len(W),tetrahedra=len(T),nodes=len(P),originalFaces=len(F),
        minimumPrismVolume=float(min(volume)),minimumTetVolume=float(min(tv)) if len(tv) else None,
        maximumPrismClosureRelative=closureError,totalFluidVolume=total,outerVolume=outer,bodyVolume=body,
        fluidVolumeClosureRelative=relative,triangular=tr,quadrilateral=qr,physicsQualified=False,newQualifiedArea=0)
