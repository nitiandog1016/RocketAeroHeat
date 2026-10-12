function testSurfaceStreamlines()
% 平面三角网格上，流动方向、连通性和距离应沿 +X 增长。
addpath(fileparts(fileparts(mfilename('fullpath'))));
nx=8; ny=3; [x,y]=ndgrid(0:nx,0:ny);
m.V=[x(:),y(:),zeros(numel(x),1)];
index=reshape(1:numel(x),size(x));
F=zeros(2*nx*ny,3); k=0;
for ix=1:nx
    for iy=1:ny
        a=index(ix,iy); b=index(ix+1,iy);
        c=index(ix+1,iy+1); d=index(ix,iy+1);
        k=k+1; F(k,:)=[a b c];
        k=k+1; F(k,:)=[a c d];
    end
end
m.F=F; nf=size(F,1);
m.localGeometry.center=(m.V(F(:,1),:)+m.V(F(:,2),:)+m.V(F(:,3),:))/3;
m.localGeometry.normal=repmat([0 0 1],nf,1);
m.boundaryLayerEdge.velocityVector=repmat([100 0 0],nf,1);
m.boundaryLayerEdge.method='解析平面均匀流';
m.stagnation.candidateFace=m.localGeometry.center(:,1)<0.5;
m=traceSurfaceStreamlines(m);
s=m.surfaceStreamlines;
assert(all(s.connectedFace));
assert(all(s.distance(s.seedFace)==0));
assert(all(abs(s.tangentDirection(:,1)-1)<1e-12));
assert(min(s.distance(m.localGeometry.center(:,1)>7))>6);
assert(max(s.distance)<12);
assert(all(s.pathCost(s.connectedFace)>=s.distance(s.connectedFace)-1e-10));
assert(any(s.pathCost(s.connectedFace)>s.distance(s.connectedFace)+0.01));
children=find(s.predecessorFace>0);
assert(all(abs(s.distance(children)-s.distance(s.predecessorFace(children)) ...
    -s.stepLength(children))<1e-10));
assert(contains(s.method,'Dijkstra'));
far=find(s.connectedFace & m.localGeometry.center(:,1)>7,1);
curve=traceSurfacePath(m,far);
assert(curve.available && size(curve.points,2)==3);
assert(abs(curve.distance-s.distance(far))<1e-10);
assert(all(diff(curve.faceIDs)~=0));
assert(all(diff(curve.points(:,1))>=-1e-12));
bad=m; bad.boundaryLayerEdge.velocityVector(:)=0;
failed=false; try, traceSurfaceStreamlines(bad); catch, failed=true; end
assert(failed);
fprintf('Surface streamline test passed: %d connected faces, maximum distance %.3f m.\n', ...
    nnz(s.connectedFace),max(s.distance));
end
