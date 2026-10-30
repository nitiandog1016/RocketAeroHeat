function mesh=buildShockFittedTetraTopology(nodes,tetrahedra)
%BUILDSHOCKFITTEDTETRATOPOLOGY 建立真正共享的三维四面体有限体积拓扑。
% 模块一：每个内部三角面只存一份，owner/neighbour保留相反有向面积。
% 模块二：拒绝负体积、重复单元、非流形面和相邻面朝向错误，不自动修补输入。
% 模块三：稀疏关联矩阵只装配±1拓扑，不建立五量Jacobian或大稠密缓存。
validateattributes(nodes,{'double'},{'real','finite','ncols',3});
validateattributes(tetrahedra,{'double'},{'real','finite','integer','ncols',4,'positive','<=',size(nodes,1)});
assert(size(nodes,1)>=4&&size(tetrahedra,1)>0&& ...
    size(unique(sort(tetrahedra,2),'rows'),1)==size(tetrahedra,1), ...
    'ShockFit:Topology','需要非空且不重复的四面体。');
a=nodes(tetrahedra(:,1),:);b=nodes(tetrahedra(:,2),:);c=nodes(tetrahedra(:,3),:);d=nodes(tetrahedra(:,4),:);
volumes=sum((b-a).*cross(c-a,d-a,2),2)/6;
assert(all(volumes>0&isfinite(volumes)),'ShockFit:CellVolume','输入单元体积须严格为正。');
Nc=size(tetrahedra,1);allFaces=[tetrahedra(:,[2 3 4]);tetrahedra(:,[1 4 3]);tetrahedra(:,[1 2 4]);tetrahedra(:,[1 3 2])];
cells=repmat((1:Nc)',4,1);[~,first,index]=unique(sort(allFaces,2),'rows');counts=accumarray(index,1);
assert(all(counts==1|counts==2),'ShockFit:NonManifold','出现非流形面，不能用于区域计算。');
faces=allFaces(first,:);owner=cells(first);neighbour=zeros(size(owner));
remaining=true(size(index));remaining(first)=false;
neighbour(index(remaining))=cells(remaining);assert(all(neighbour~=owner),'ShockFit:Topology','同一单元重复使用同一面。');
Nf=size(faces,1);interior=neighbour>0;
incidence=sparse([owner;neighbour(interior)],[(1:Nf)';find(interior)], ...
    [ones(Nf,1);-ones(nnz(interior),1)],Nc,Nf);
assert(all(full(sum(abs(incidence),2))==4),'ShockFit:Topology','四面体没有恰好四个面。');
mesh=struct('tetrahedra',tetrahedra,'faces',faces,'owner',owner,'neighbour',neighbour, ...
    'incidence',incidence,'nodeCount',size(nodes,1),'cellCount',Nc,'faceCount',Nf,'boundaryFaces',find(~interior));
geo=shockFittedALEGeometry(mesh,nodes,nodes,1);
assert(max(geo.areaClosureRelative)<2e-12,'ShockFit:Topology','有向面没有围成闭合单元。');
end
