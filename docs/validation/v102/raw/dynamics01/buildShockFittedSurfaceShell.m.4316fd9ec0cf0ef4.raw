function shell=buildShockFittedSurfaceShell(V,F,originalFaceIDs,gapM,layers)
%BUILDSHOCKFITTEDSURFACESHELL 从真实原面元生成单一相连的三维局部体网格。
% 模块一：壁面节点和三角面逐字保持原几何；外表面为共享顶点法向偏移候选。
% 模块二：所有棱柱按全局节点编号一致剖分，内层/邻面共享，不复制轴对称剖面。
% 模块三：偏移外壳尚不是求得的激波；保存原面号，不能据此授予外缘/热流资格。
validateattributes(V,{'double'},{'real','finite','ncols',3});
validateattributes(F,{'double'},{'real','finite','integer','positive','ncols',3,'<=',size(V,1)});
validateattributes(originalFaceIDs,{'double'},{'real','finite','integer','positive','column','<=',size(F,1)});
validateattributes(gapM,{'double'},{'real','finite','positive','scalar'});
validateattributes(layers,{'double'},{'real','finite','integer','positive','scalar'});
assert(numel(unique(originalFaceIDs))==numel(originalFaceIDs),'ShockFit:OriginalFaces','原面号重复。');
patch=F(originalFaceIDs,:);[ids,~,local]=unique(patch(:));patch=reshape(local,size(patch));body=V(ids,:);
vectors=cross(body(patch(:,2),:)-body(patch(:,1),:),body(patch(:,3),:)-body(patch(:,1),:),2);
normals=zeros(size(body));
for j=1:3,for k=1:3,normals(:,k)=normals(:,k)+accumarray(patch(:,j),vectors(:,k),[numel(ids) 1]);end,end
assert(all(vecnorm(normals,2,2)>0),'ShockFit:ShellNormals','共享顶点法向退化。');
normals=normals./vecnorm(normals,2,2);N=size(body,1);nodes=zeros((layers+1)*N,3);
for layer=0:layers,nodes(layer*N+(1:N),:)=body+(layer/layers)*gapM*normals;end
sorted=sort(patch,2);tetra=zeros(layers*3*size(patch,1),4);sourceFace=zeros(size(tetra,1),1);q=0;
for layer=0:layers-1
    a=sorted(:,1)+layer*N;b=sorted(:,2)+layer*N;c=sorted(:,3)+layer*N;
    group=[a b c c+N;a b b+N c+N;a a+N b+N c+N];
    index=q+(1:size(group,1));tetra(index,:)=group;sourceFace(index)=repmat(originalFaceIDs,3,1);q=q+size(group,1);
end
a=nodes(tetra(:,1),:);b=nodes(tetra(:,2),:);c=nodes(tetra(:,3),:);d=nodes(tetra(:,4),:);
negative=sum((b-a).*cross(c-a,d-a,2),2)<0;tetra(negative,[3 4])=tetra(negative,[4 3]);
mesh=buildShockFittedTetraTopology(nodes,tetra);boundaryNodes=mesh.faces(mesh.boundaryFaces,:);
wall=all(boundaryNodes<=N,2);outer=all(boundaryNodes>layers*N,2);
assert(nnz(wall)==size(patch,1)&&nnz(outer)==size(patch,1),'ShockFit:ShellTopology','真实壁面或外壳面数变化。');
shell=struct('nodes',nodes,'mesh',mesh,'originalVertexIDs',ids,'bodyFaces',patch,'originalFaceIDs',originalFaceIDs, ...
    'cellOriginalFaceIDs',sourceFace,'vertexNormals',normals,'gapM',gapM,'layers',layers, ...
    'wallBoundaryMask',wall,'outerBoundaryMask',outer,'sideBoundaryMask',~wall&~outer, ...
    'bodyVerticesExact',isequal(nodes(1:N,:),V(ids,:)),'shockGeometrySolved',false, ...
    'pressureOnBodySolved',false,'heatFluxSolved',false,'wallTemperatureSolved',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0);
end
