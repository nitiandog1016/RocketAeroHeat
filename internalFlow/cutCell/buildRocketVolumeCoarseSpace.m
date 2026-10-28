function space=buildRocketVolumeCoarseSpace(cut,options)
%BUILDROCKETVOLUMECOARSESPACE 为原面法向层域构造体积加权、连通聚合的代数粗空间。
% 模块一：坐标分箱和法向层分段仅定义预条件聚合，不生成或替换物理网格。
% 模块二：只沿正面积共享面连接同箱单元，分开不连通分量，不丢弃小控制体。
% 模块三：P为分片常数延拓，R为真实V加权平均；核对R*P=I与总体积。
if nargin<2,options=struct;end
nf=size(cut.inputGeometry.F,1);n=numel(cut.cells.volume);nl=n/nf;
assert(strcmp(cut.kind,'original-triangle-body-fitted-three-dimensional-shell')&&nl==floor(nl),'RocketCoarse:Geometry','当前粗空间要求原面法向层域。');
defaults=struct('cartesianBins',[2 2 2],'layerBreaks',floor(nl/2),'maximumCoarseUnknowns',200);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'RocketCoarse:Options','未知粗空间选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.cartesianBins,{'numeric'},{'real','finite','integer','positive','numel',3});
validateattributes(options.maximumCoarseUnknowns,{'numeric'},{'real','finite','integer','positive','scalar'});
breaks=options.layerBreaks(:)';assert(isnumeric(breaks)&&isreal(breaks)&&all(isfinite(breaks))&&all(breaks==floor(breaks)&breaks>=1&breaks<nl)&&all(diff(breaks)>0), ...
    'RocketCoarse:Layers','层分段须为严格递增的内部层编号。');
v=cut.cells.volume(:);c=cut.cells.centroid;assert(all(v>0&isfinite(v))&&all(isfinite(c),'all'));
bin=zeros(n,3);thresholds=cell(1,3);
for d=1:3
    sorted=sort(c(:,d));count=options.cartesianBins(d);indices=max(1,min(n,round((1:count-1)*n/count)));
    thresholds{d}=unique(sorted(indices));
    for t=thresholds{d}',bin(:,d)=bin(:,d)+double(c(:,d)>t);end
end
layer=floor((0:n-1)'/nf)+1;radial=zeros(n,1);for t=breaks,radial=radial+double(layer>t);end
[~,~,bucket]=unique([bin radial],'rows');f=cut.openFaces;
inside=f.area>0&f.cellMinus>0&f.cellPlus>0;m=double(f.cellMinus(inside));p=double(f.cellPlus(inside));same=bucket(m)==bucket(p);
% 对称正邻接加对角后，dmperm的不可约块对应每个箱内的连通分量。
A=sparse([m(same);p(same);(1:n)'],[p(same);m(same);(1:n)'],1,n,n);
[order,~,boundaries,~]=dmperm(A);groups=numel(boundaries)-1;map=zeros(n,1);
for g=1:groups,map(order(boundaries(g):boundaries(g+1)-1))=g;end
assert(5*groups<=options.maximumCoarseUnknowns,'RocketCoarse:Budget','连通粗空间有%d未知量，超过预先规定的%d；不截掉分量。',5*groups,options.maximumCoarseUnknowns);
groupVolume=accumarray(map,v,[groups 1]);P=sparse((1:n)',map,1,n,groups);R=sparse(map,(1:n)',v./groupVolume(map),groups,n);
identityError=norm(R*P-speye(groups),Inf);volumeError=abs(sum(groupVolume)-sum(v))/sum(v);
assert(identityError<1e-12&&volumeError<1e-12,'RocketCoarse:Transfer','体积加权传递或体积账不一致。');
space=struct('kind','volume-weighted-connected-original-CV-aggregation-for-preconditioning','geometry',cut, ...
    'cellToGroup',map,'fineVolume',v,'groupVolume',groupVolume,'groupCentroid',R*c,'prolongation',P,'restriction',R, ...
    'groupCount',groups,'coarseUnknowns',5*groups,'coordinateThresholds',{thresholds},'options',options, ...
    'restrictionProlongationIdentityError',identityError,'totalVolumeRelativeDifference',volumeError, ...
    'physicalMeshChanged',false,'newProductionSolvedAreaM2',0);
end
