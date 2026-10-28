function [fv,domain]=buildRocketAxisymmetricCollarDomain(surface,truth,domainOptions,collarOptions)
%BUILDROCKETAXISYMMETRICCOLLARDOMAIN 完整箭体的分段等法向偏移近壁层研究网格。
% 模块一：复用原构造器严格核对真实完整母线、源面、两轴边和远场，不改表面网格。
% 模块二：每个节点求相邻墙线的等距交点，鼻尖/底心的偏移节点仍严格位于轴上。
% 模块三：前三层按真实墙法向距离偏移，再以显式几何增长直线接回同一完整远场。
% 模块四：重新计算全部FV度量、WLS和墙距，拒绝折叠/越界，不裁网格或流动状态。
% 这是未接入生产流场/源面场量投影的几何候选，不建立稳态、y+或物理求解资格。
if nargin<1||isempty(surface)
    [surface,truth]=buildRocketTestBody(struct('protrusionEnabled',false));
elseif nargin<2||isempty(truth)
    error('RocketCollar:Source','显式表面必须附带原生完整几何真值。');
end
if nargin<3,domainOptions=struct('radialLayers',24);end
if nargin<4,collarOptions=struct;end
defaults=struct('normalLayers',3,'firstNormalSpacing',2e-5,'normalGrowth',1.5);
assert(isstruct(collarOptions)&&isscalar(collarOptions)&& ...
    all(ismember(fieldnames(collarOptions),fieldnames(defaults))), ...
    'RocketCollar:Options','只接受近壁层数、首法向偏移及法向增长率。');
for key=fieldnames(defaults)'
    name=key{1};if ~isfield(collarOptions,name),collarOptions.(name)=defaults.(name);end
end
validateattributes(collarOptions.normalLayers,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(collarOptions.firstNormalSpacing,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(collarOptions.normalGrowth,{'numeric'},{'real','finite','scalar','>=',1});
[legacy,domain]=buildRocketAxisymmetricDomain(surface,truth,domainOptions);
n=domain.options.radialLayers;nc=collarOptions.normalLayers;
assert(nc<=n-2,'RocketCollar:Layers','法向近壁层外至少保留两层连接远场。');
wall=domain.wallPolylineXR;outer=domain.outerBoundaryPolylineXR;edge=diff(wall);
normal=[edge(:,2),-edge(:,1)]./vecnorm(edge,2,2);nn=size(wall,1);direction=zeros(nn,2);
assert(normal(1,1)>0&&normal(end,1)<0,'RocketCollar:Axis','真实鼻尖/底心偏移必须沿各自的外向轴。');
direction(1,1)=1/normal(1,1);direction(end,1)=1/normal(end,1);

%% 一、真实相邻墙线等距交点；共线节点使用角平分法，反向/退化角拒绝
for k=2:nn-1
    a=normal(k-1,:);b=normal(k,:);den=1+dot(a,b);
    assert(isfinite(den)&&den>0,'RocketCollar:Corner','相邻墙线接近反向，不能建立有限等距交点。');
    direction(k,:)=(a+b)/den;
end
left=sum(direction(1:end-1,:).*normal,2);right=sum(direction(2:end,:).*normal,2);
budget=512*eps*(1+sum(abs(direction(1:end-1,:).*normal),2)+sum(abs(direction(2:end,:).*normal),2));
assert(all(isfinite(direction),'all')&&all(abs(left-1)<=budget&abs(right-1)<=budget), ...
    'RocketCollar:Corner','偏移方向未满足两相邻真实墙线的单位法向约束。');
offset=[0,cumsum(collarOptions.firstNormalSpacing*collarOptions.normalGrowth.^(0:nc-1))];
assert(all(isfinite(offset))&&all(diff(offset)>0),'RocketCollar:Spacing','法向层偏移须有限严格递增。');
X=zeros(nn,n+1);R=X;
X(:,1:nc+1)=wall(:,1)+direction(:,1)*offset;
R(:,1:nc+1)=wall(:,2)+direction(:,2)*offset;
start=[X(:,nc+1),R(:,nc+1)];connection=outer-start;span=vecnorm(connection,2,2);
assert(all(isfinite(start),'all')&&all(start(:,2)>=0)&& ...
    all(vecnorm(start-domain.kernel,2,2)<domain.options.outerRadius), ...
    'RocketCollar:Extent','法向层越过轴或远场；不得裁点或静默缩减给定偏移。');
remaining=n-nc;first=domain.options.firstRadialSpacing;
assert(all(isfinite(span)&span>=remaining*first), ...
    'RocketCollar:Extent','给定远场不足以容纳显式外部连接层。');
growth=zeros(nn,1);

%% 二、从近壁层外缘接回原远场；该首连接长度不等于墙法向距离
for k=1:nn
    growth(k)=growthRatio(first,span(k),remaining);
    fraction=cumsum(first*growth(k).^(0:remaining-1))/span(k);fraction(end)=1;
    assert(all(isfinite(fraction))&&all(diff([0 fraction])>0)&&all(fraction<=1), ...
        'RocketCollar:Spacing','外部连接层必须严格向给定远场推进。');
    X(k,nc+2:end)=start(k,1)+connection(k,1)*fraction;
    R(k,nc+2:end)=start(k,2)+connection(k,2)*fraction;
end
X(:,1)=wall(:,1);R(:,1)=wall(:,2);X(:,end)=outer(:,1);R(:,end)=outer(:,2);
R(1,:)=0;R(end,:)=0;
fv=buildAxisymmetricFVGeometry(X,R);
assert(isequal(fv.faces.owner,legacy.faces.owner)&&isequal(fv.faces.neighbour,legacy.faces.neighbour)&& ...
    isequal(fv.faces.boundarySide,legacy.faces.boundarySide),'RocketCollar:Topology','新网格必须保持同一完整C域连接关系。');
ids=domain.wallBoundaryFaceID;axes=[domain.axisBoundaryFaceID{1};domain.axisBoundaryFaceID{2}];
assert(isequal(fv.faces.centroid(ids,:),legacy.faces.centroid(ids,:))&& ...
    isequal(fv.faces.areaPerRadian(ids),legacy.faces.areaPerRadian(ids))&& ...
    all(fv.faces.centroid(axes,2)==0)&all(fv.faces.areaPerRadian(axes)==0), ...
    'RocketCollar:Boundary','不能移动真实墙、改变墙面积或把轴移出r=0。');
volume=sum(fv.volumePerRadian);expected=domain.outerPolygonVolumePerRadianM3-domain.bodyVolumePerRadianM3;
assert(abs(volume-expected)<=2e-12*(abs(volume)+abs(expected)), ...
    'RocketCollar:Volume','新网格完整r权体积与同一外域减实体不符。');
legacySkew=skew(legacy,ids);newSkew=skew(fv,ids);

%% 三、所有依赖新网格的度量重新绑定；源面候选保留，尚不宣称场量投影完成
domain.kind='complete-native-rocket-normal-collar-exterior-domain-v1';
domain.meshFactory='buildRocketAxisymmetricCollarDomain';
domain.radialGrowthRatio=growth;
domain.radialGrowthRatioMeaning='outer straight connection only; inner layers use separate normalGrowth';
domain.wallOwnerCell=fv.faces.owner(ids);domain.wallNormalDistanceM=fv.faces.normalDistance(ids);
domain.axisFaceCentroidXR=fv.faces.centroid(axes,:);
domain.axisFaceOwnerXOffsetM=fv.faces.centroid(axes,1)-fv.volumeCentroid(fv.faces.owner(axes),1);
domain.fluidVolumePerRadianM3=volume;domain.metricClosureRelative=fv.metricClosureRelative;
domain.collar=struct('options',collarOptions,'vertexOffsetDirectionXR',direction, ...
    'wallOutwardNormalXR',normal,'normalOffsetM',offset,'outerConnectionFirstLengthM',first, ...
    'perpendicularVertexOffsetNotCellCentroidDistance',true,'firstNormalSpacingIsOwnerWallDistance',false, ...
    'legacyFirstCellNormalDistanceM',legacy.faces.normalDistance(ids),'legacyFirstCellSkewRatio',legacySkew, ...
    'firstCellSkewRatio',newSkew,'completeWallAndOuterBoundaryPreserved',true, ...
    'surfaceRemesherChanged',false,'sourceFieldProjectionCompleted',false,'flowAdvancePerformed',false, ...
    'nearWallResolutionValidated',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
end

function value=skew(mesh,ids)
% 模块用途：真实r权面心至体质心的切向偏移/法向墙距，用于几何对照而非y+验收。
t=[-mesh.faces.normal(ids,2),mesh.faces.normal(ids,1)];
value=abs(sum(mesh.faces.delta(ids,:).*t,2))./mesh.faces.normalDistance(ids);
end
function q=growthRatio(first,span,count)
% 模块用途：给定真实路径总长及首长度，求外部连接层的非减几何增长率。
if abs(span-first*count)<=128*eps*max(span,first*count),q=1;return;end
lo=1;hi=2;
while first*sum(hi.^(0:count-1))<span
    hi=2*hi;assert(isfinite(hi),'RocketCollar:Spacing','外部层增长率不可表示。');
end
for k=1:80
    mid=(lo+hi)/2;
    if first*sum(mid.^(0:count-1))<span,lo=mid;else,hi=mid;end
end
q=(lo+hi)/2;
end
