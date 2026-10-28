function space=buildRocketTargetedCoarseSpace(parent,options)
%BUILDROCKETTARGETEDCOARSESPACE 只在鼻端近壁细分旧连通粗组，不改变原控制体或空间方程。
% 模块一：按原墙面轴向位置选目标面，在指定法向层内按CV质心X中位数再分组。
% 模块二：保留旧父组边界，沿正面积共享面分开不连通部分；小分量不删除。
% 模块三：真实体积R/P与父子嵌套逐项核对，退化的单点/同坐标组保持原样。
if nargin<2,options=struct;end
cut=parent.geometry;nf=size(cut.inputGeometry.F,1);n=numel(cut.cells.volume);nl=n/nf;
defaults=struct('axialFractionThreshold',.9,'targetLayers',1:min(4,nl),'maximumCoarseUnknowns',500);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'TargetCoarse:Options','未知局部细分选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.axialFractionThreshold,{'numeric'},{'real','finite','scalar','>',0,'<',1});
validateattributes(options.maximumCoarseUnknowns,{'numeric'},{'real','finite','scalar','integer','positive'});
layers=options.targetLayers(:)';assert(isnumeric(layers)&&all(isfinite(layers))&&all(layers==floor(layers)&layers>=1&layers<=nl)&&all(diff(layers)>0),'TargetCoarse:Layers','目标层须为严格递增的合法编号。');
V=cut.inputGeometry.V;F=cut.inputGeometry.F;wallX=mean(reshape(V(F(:),1),nf,3),2);
fraction=(wallX-min(V(:,1)))/(max(V(:,1))-min(V(:,1)));layer=floor((0:n-1)'/nf)+1;
target=repmat(fraction>=options.axialFractionThreshold,nl,1)&ismember(layer,layers');
bucket=zeros(n,1);thresholds=nan(parent.groupCount,1);
for g=1:parent.groupCount
    ids=find(parent.cellToGroup==g&target);if isempty(ids),continue;end
    x=cut.cells.centroid(ids,1);thresholds(g)=median(x);bucket(ids)=1+double(x>thresholds(g));
end
f=cut.openFaces;inside=f.area>0&f.cellMinus>0&f.cellPlus>0;m=double(f.cellMinus(inside));p=double(f.cellPlus(inside));
same=parent.cellToGroup(m)==parent.cellToGroup(p)&bucket(m)==bucket(p);
A=sparse([m(same);p(same);(1:n)'],[p(same);m(same);(1:n)'],1,n,n);[order,~,boundaries,~]=dmperm(A);
ng=numel(boundaries)-1;assert(5*ng<=options.maximumCoarseUnknowns,'TargetCoarse:Budget','局部连通细分%d未知量超过%d，不删除分量。',5*ng,options.maximumCoarseUnknowns);
map=zeros(n,1);for g=1:ng,map(order(boundaries(g):boundaries(g+1)-1))=g;end
v=cut.cells.volume(:);vol=accumarray(map,v,[ng 1]);P=sparse((1:n)',map,1,n,ng);R=sparse(map,(1:n)',v./vol(map),ng,n);
transfer=R*parent.prolongation;identity=norm(R*P-speye(ng),Inf);nesting=norm(P*transfer-parent.prolongation,Inf);
assert(identity<1e-12&&nesting<1e-12&&abs(sum(vol)-sum(v))/sum(v)<1e-12,'TargetCoarse:Transfer','体积传递或父空间嵌套失败。');
space=parent;space.cellToGroup=map;space.fineVolume=v;space.groupVolume=vol;space.groupCentroid=R*cut.cells.centroid;
space.prolongation=P;space.restriction=R;space.groupCount=ng;space.coarseUnknowns=5*ng;
space.restrictionProlongationIdentityError=identity;space.totalVolumeRelativeDifference=abs(sum(vol)-sum(v))/sum(v);
space.targetRefinement=struct('options',options,'targetCellCount',sum(target),'groupXThresholds',thresholds, ...
    'parentGroupCount',parent.groupCount,'parentCellToGroup',parent.cellToGroup,'parentTransfer',transfer,'parentNestingInfinityError',nesting);
end
