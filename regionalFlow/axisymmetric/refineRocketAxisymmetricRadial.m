function [mesh,domain,audit]=refineRocketAxisymmetricRadial(parent,parentDomain)
%REFINEROCKETAXISYMMETRICRADIAL 逐原控制体径向二分，保留真实墙和外边界。
% 模块一：每相邻径向节点列取中点，不移动父节点或重新拟合完整母线。
% 模块二：完整重建r权体积、两种质心、共享面、WLS和墙/轴/外边界引用。
% 模块三：逐父格检查两个子格的体积及r权一阶矩相加；不是流场保守重映射。
assert(isequal(parent.shape,[size(parent.X,1)-1,size(parent.X,2)-1])&& ...
    isequal([parent.X(:,1),parent.R(:,1)],parentDomain.wallPolylineXR), ...
    'RocketRadialRefine:Parent','需要完整且匹配的贴体父网格。');
n=parent.shape(2);X=zeros(size(parent.X,1),2*n+1);R=X;
X(:,1:2:end)=parent.X;R(:,1:2:end)=parent.R;
X(:,2:2:end)=.5*(parent.X(:,1:end-1)+parent.X(:,2:end));
R(:,2:2:end)=.5*(parent.R(:,1:end-1)+parent.R(:,2:end));
mesh=buildAxisymmetricFVGeometry(X,R);ni=parent.shape(1);
id=reshape(1:mesh.cellCount,ni,2*n);a=id(:,1:2:end);b=id(:,2:2:end);a=a(:);b=b(:);
v=mesh.volumePerRadian(a)+mesh.volumePerRadian(b);mom=mesh.volumePerRadian(a).*mesh.volumeCentroid(a,:)+ ...
    mesh.volumePerRadian(b).*mesh.volumeCentroid(b,:);oldmom=parent.volumePerRadian.*parent.volumeCentroid;
relative=max(abs(v-parent.volumePerRadian)./parent.volumePerRadian);
momentError=max(abs(mom-oldmom)./(parent.volumePerRadian.*max(1,abs(parent.volumeCentroid))),[],'all');
assert(relative<3e-10&&momentError<3e-10,'RocketRadialRefine:Moments','子格真实体积/一阶矩未还原父格。');
domain=parentDomain;wall=find(mesh.faces.boundarySide==3);ordered=zeros(ni,1);
for i=1:ni,match=wall(mesh.faces.owner(wall)==i);assert(isscalar(match));ordered(i)=match;end
domain.wallBoundaryFaceID=ordered;domain.wallOwnerCell=mesh.faces.owner(ordered);
domain.wallNormalDistanceM=mesh.faces.normalDistance(ordered);
domain.axisBoundaryFaceID={find(mesh.faces.boundarySide==1),find(mesh.faces.boundarySide==2)};
axis=[domain.axisBoundaryFaceID{1};domain.axisBoundaryFaceID{2}];
domain.axisFaceCentroidXR=mesh.faces.centroid(axis,:);
domain.axisFaceOwnerXOffsetM=mesh.faces.centroid(axis,1)-mesh.volumeCentroid(mesh.faces.owner(axis),1);
domain.outerBoundaryFaceID=find(mesh.faces.boundarySide==4);
domain.options.radialLayers=2*n;domain.options.firstRadialSpacing=parentDomain.options.firstRadialSpacing/2;
domain.kind='complete-native-rocket-nested-radial-bisection-domain-v1';domain.meshFactory='refineRocketAxisymmetricRadial';
domain.fluidVolumePerRadianM3=sum(mesh.volumePerRadian);domain.metricClosureRelative=mesh.metricClosureRelative;
domain.radialGrowthRatioMeaning='parent outer-layer growth; child thickness sequence repeats half of each parent thickness, not a single geometric series';
if isfield(domain,'collar')
    old=parentDomain.collar;off=old.normalOffsetM;offset=zeros(1,2*numel(off)-1);offset(1:2:end)=off;offset(2:2:end)=.5*(off(1:end-1)+off(2:end));
    domain.collar.normalOffsetM=offset;domain.collar.options.normalLayers=2*old.options.normalLayers;
    domain.collar.options.firstNormalSpacing=old.options.firstNormalSpacing/2;
    domain.collar.options.radialRefinementFactor=2;
    domain.collar.options.normalGrowthMeaning='parent-layer growth; child consecutive thickness ratios alternate 1 and parent growth';
    domain.collar.parentRequestedOptions=old.options;
    domain.collar.outerConnectionFirstLengthM=old.outerConnectionFirstLengthM/2;
    t=[-mesh.faces.normal(ordered,2),mesh.faces.normal(ordered,1)];
    domain.collar.firstCellSkewRatio=abs(sum(mesh.faces.delta(ordered,:).*t,2))./domain.wallNormalDistanceM;
    domain.collar.legacyReferenceIsParentOnly=true;
end
assert(isequal(mesh.X(:,1:2:end),parent.X)&&isequal(mesh.R(:,1:2:end),parent.R)&& ...
    isequal(mesh.faces.centroid(ordered,:),parent.faces.centroid(parentDomain.wallBoundaryFaceID,:))&& ...
    isequal(mesh.faces.areaPerRadian(ordered),parent.faces.areaPerRadian(parentDomain.wallBoundaryFaceID)), ...
    'RocketRadialRefine:Boundary','二分不能移动父节点、墙或墙面积。');
audit=struct('kind','exact-parent-node-preserving-radial-bisection-v1','parentCellCount',parent.cellCount, ...
    'childCellCount',mesh.cellCount,'parentToChildren',[a b], ...
    'maximumParentVolumeRelativeError',relative,'maximumScaledFirstMomentErrorM',momentError, ...
    'wallAndOuterBoundaryPreserved',true,'streamwiseRefinementPerformed',false, ...
    'conservativeFlowRemappingPerformed',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
domain.radialRefinement=audit;
end
