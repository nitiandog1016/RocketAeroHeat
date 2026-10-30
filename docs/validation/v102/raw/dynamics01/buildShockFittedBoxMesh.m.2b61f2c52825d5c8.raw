function box=buildShockFittedBoxMesh(counts,lengths)
%BUILDSHOCKFITTEDBOXMESH 生成三维共享四面体盒域，供非均匀波传播和坐标变换核验。
% 模块一：每个六面体一致剖分为6个四面体，共享面仅一份。
% 模块二：输出各边界坐标方向，不将此规则网格当作真实箭体几何。
validateattributes(counts,{'double'},{'real','finite','integer','positive','size',[1 3]});
validateattributes(lengths,{'double'},{'real','finite','positive','size',[1 3]});
[x,y,z]=ndgrid(linspace(0,lengths(1),counts(1)+1),linspace(0,lengths(2),counts(2)+1),linspace(0,lengths(3),counts(3)+1));
X=[x(:) y(:) z(:)];ids=reshape(1:numel(x),size(x));T=zeros(6*prod(counts),4);q=0;
for k=1:counts(3),for j=1:counts(2),for i=1:counts(1)
    a=ids(i,j,k);b=ids(i+1,j,k);c=ids(i+1,j+1,k);d=ids(i,j+1,k);
    e=ids(i,j,k+1);f=ids(i+1,j,k+1);g=ids(i+1,j+1,k+1);h=ids(i,j+1,k+1);
    T(q+(1:6),:)=[a b c g;a c d g;a d h g;a h e g;a e f g;a f b g];q=q+6;
end,end,end
v=sum((X(T(:,2),:)-X(T(:,1),:)).*cross(X(T(:,3),:)-X(T(:,1),:),X(T(:,4),:)-X(T(:,1),:),2),2);
negative=v<0;T(negative,[3 4])=T(negative,[4 3]);mesh=buildShockFittedTetraTopology(X,T);
centers=(X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:)+X(T(:,4),:))/4;
faces=mesh.faces(mesh.boundaryFaces,:);fc=(X(faces(:,1),:)+X(faces(:,2),:)+X(faces(:,3),:))/3;
left=all(reshape(X(faces(:),1),size(faces))==0,2);
right=all(reshape(X(faces(:),1),size(faces))==lengths(1),2);
box=struct('nodes',X,'mesh',mesh,'cellCentroidM',centers,'boundaryCentroidM',fc, ...
    'leftBoundaryMask',left,'rightBoundaryMask',right,'counts',counts,'lengthsM',lengths);
end
