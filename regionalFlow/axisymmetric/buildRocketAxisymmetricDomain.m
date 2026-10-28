function [fv,domain]=buildRocketAxisymmetricDomain(surface,truth,options)
%BUILDROCKETAXISYMMETRICDOMAIN 为完整无突起箭体建立轴对称外流C域。
% 模块一：核对原生闭合表面及精确圆半径，保留鼻尖、接缝和真实底面。
% 模块二：从实体内核点按严格递增极角连接同射线远场，生成正凸四边形。
% 模块三：每条射线使用显式首径向间距及独立几何增长率，不冒称法向/y+合格。
% 模块四：保存真实墙段、原faceBand/三角面候选、轴边和完整面积/体积账。
% 本函数只构建真实几何；偏斜轴邻近须启用独立镜像模板，不能删旧AxisSkew护栏。
% 不支持非轴对称突起，不输入外部CFD，不建立任何新增物理解面积资格。
if nargin<1||isempty(surface)
    [surface,truth]=buildRocketTestBody(struct('protrusionEnabled',false));
elseif nargin<2||isempty(truth)
    error('RocketDomain:Source','显式表面必须附带原生几何真值。');
end
if nargin<3,options=struct;end
assert(isstruct(options)&&isscalar(options),'RocketDomain:Options','options须为标量结构体。');
defaults=struct('kernelX',[],'outerRadius',[],'radialLayers',12,'firstRadialSpacing',2e-4, ...
    'maxWallAngleStepDeg',5,'baseRadialSegments',8);
assert(all(ismember(fieldnames(options),fieldnames(defaults))),'RocketDomain:Options','存在未知域参数。');
for key=fieldnames(defaults)',k=key{1};if ~isfield(options,k),options.(k)=defaults.(k);end,end
validateattributes(options.radialLayers,{'numeric'},{'real','finite','scalar','integer','>=',2,'<=',512});
validateattributes(options.baseRadialSegments,{'numeric'},{'real','finite','scalar','integer','>=',1,'<=',512});
validateattributes(options.firstRadialSpacing,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(options.maxWallAngleStepDeg,{'numeric'},{'real','finite','scalar','>',0,'<=',45});

%% 一、原生完整表面合同；不从多边形等面积半径替代真实母线
assert(isstruct(surface)&&isscalar(surface)&&all(isfield(surface,{'V','F','param'}))&& ...
    isstruct(truth)&&isscalar(truth)&&all(isfield(truth,{'kind','features','nodeBaselineRadius', ...
    'nodeRadialOffset','faceBaseRegion','totalLength','cylinderX','baseRadius','exactBaselineArea','exactBaselineVolume'}))&& ...
    strcmp(truth.kind,'parametricRocketRadialGraph'),'RocketDomain:Source','需要buildRocketTestBody原生完整表面及真值。');
assert(isempty(truth.features)&&all(truth.nodeRadialOffset==0),'RocketDomain:Nonaxisymmetric','轴对称域明确拒绝突起/偏心截面。');
V=surface.V;F=surface.F;p=surface.param;
assert(isnumeric(V)&&isreal(V)&&size(V,2)==3&&all(isfinite(V),'all')&& ...
    isnumeric(F)&&isreal(F)&&size(F,2)==3&&all(isfinite(F)&F==fix(F)&F>=1&F<=size(V,1),'all')&& ...
    all(isfield(p,{'sectionX','sectionNodeIDs','noseCapNodeID','aftCapNodeID','faceBandID'})), ...
    'RocketDomain:Source','表面节点/面元或原生截面合同不完整。');
assert(numel(truth.nodeBaselineRadius)==size(V,1)&&numel(truth.nodeRadialOffset)==size(V,1)&& ...
    numel(truth.faceBaseRegion)==size(F,1)&&numel(p.faceBandID)==size(F,1), ...
    'RocketDomain:Source','几何真值必须逐节点/逐面绑定。');
scale=max(1,max(abs(V),[],'all'));tol=512*eps(scale);
sx=p.sectionX(:);rings=p.sectionNodeIDs(:);nr=numel(rings);
assert(nr>=3&&numel(sx)==nr&&all(isfinite(sx))&&all(diff(sx)>0)&&sx(1)==0, ...
    'RocketDomain:Source','完整母线必须从真实底面起按X递增。');
radius=zeros(nr,1);representative=zeros(nr,1);allRingIDs=[];
for i=1:nr
    ids=rings{i}(:);assert(numel(ids)>=3&&all(ids==fix(ids)&ids>=1&ids<=size(V,1))&&numel(unique(ids))==numel(ids), ...
        'RocketDomain:Source','截面节点编号无效或不闭合。');
    rr=hypot(V(ids,2),V(ids,3));radius(i)=truth.nodeBaselineRadius(ids(1));representative(i)=ids(1);
    assert(radius(i)>0&&all(abs(V(ids,1)-sx(i))<=tol)&&all(abs(rr-radius(i))<=tol)&& ...
        all(abs(truth.nodeBaselineRadius(ids)-radius(i))<=tol), ...
        'RocketDomain:Nonaxisymmetric','截面必须为轴心在原点的同一真实圆半径。');
    allRingIDs=[allRingIDs;ids]; %#ok<AGROW>
end
tip=p.noseCapNodeID;cap=p.aftCapNodeID;
assert(isscalar(tip)&&isscalar(cap)&&tip~=cap&&all([tip cap]>=1&[tip cap]<=size(V,1)&[tip cap]==fix([tip cap]))&& ...
    isequal(sort([allRingIDs;tip;cap]),(1:size(V,1))')&&numel(unique(allRingIDs))==numel(allRingIDs)&& ...
    isequal(V(tip,2:3),[0 0])&&isequal(V(cap,:),[0 0 0])&& ...
    V(tip,1)>sx(end)&&abs(V(tip,1)-truth.totalLength)<=tol&&abs(radius(1)-truth.baseRadius)<=tol, ...
    'RocketDomain:Source','鼻尖/底面中心必须是原生唯一节点，不能补虚盖或漏鼻尖。');
edges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];[~,~,eg]=unique(sort(edges,2),'rows');
counts=accumarray(eg,1);directions=accumarray(eg,sign(edges(:,2)-edges(:,1)));
adj=sparse([F(:,1);F(:,2);F(:,3)],[F(:,2);F(:,3);F(:,1)],true,size(V,1),size(V,1));
assert(all(counts==2)&&all(directions==0)&&max(conncomp(graph(adj+adj')))==1&& ...
    size(V,1)-numel(counts)+size(F,1)==2,'RocketDomain:Topology','原模型须为单一闭合一致定向实体。');
crosses=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
sourceAreas=.5*vecnorm(crosses,2,2);assert(all(sourceAreas>0), ...
    'RocketDomain:Topology','原模型含退化三角形。');
sourceVolume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
assert(sourceVolume>0,'RocketDomain:Topology','原模型朝向错误。');
original=[V(tip,1),0;sx(end:-1:1),radius(end:-1:1);0,0];
originalVertexIDs=[tip;representative(end:-1:1);cap];
originalBands=[(nr:-1:1)';-1];
assert(isequal(sort(unique(p.faceBandID(:))),[-1;(1:nr)']), ...
    'RocketDomain:Source','源面band缺失或存在未知额外面。');
originalRegions=cell(nr+1,1);
for i=1:nr+1
    ids=find(p.faceBandID(:)==originalBands(i));labels=unique(truth.faceBaseRegion(ids));
    assert(isscalar(labels),'RocketDomain:Source','一个原生母线段不能混用区域标签。');originalRegions(i)=labels;
    band=originalBands(i);triangles=F(ids,:);
    if band==nr
        valid=ismember(triangles,[rings{nr}(:);tip]);hasEnds=sum(triangles==tip,2)==1;
        expectedFaces=numel(rings{nr});
    elseif band==-1
        valid=ismember(triangles,[rings{1}(:);cap]);hasEnds=sum(triangles==cap,2)==1;
        expectedFaces=numel(rings{1});
    else
        valid=ismember(triangles,[rings{band}(:);rings{band+1}(:)]);
        hasEnds=any(ismember(triangles,rings{band}(:)),2)&any(ismember(triangles,rings{band+1}(:)),2);
        expectedFaces=numel(rings{band})+numel(rings{band+1});
    end
    assert(numel(ids)==expectedFaces&&all(valid,'all')&&all(hasEnds), ...
        'RocketDomain:Source','faceBand标签与实际源环/鼻尖/底心三角形不一致。');
end
assert(strcmp(originalRegions{end},'base'),'RocketDomain:Source','最后一段必须是真实底面。');
regionSequence=originalRegions([true;~strcmp(originalRegions(2:end),originalRegions(1:end-1))]);
assert(isequal(regionSequence,{'sphericalNose';'cone';'cylinder';'contraction';'afterbody';'base'}), ...
    'RocketDomain:Source','完整原生箭体的区域次序不符合真实母线。');

%% 二、严格星形内核与真实壁面分段；相邻射线夹角在(0,pi)内
if isempty(options.kernelX),options.kernelX=mean(truth.cylinderX);end
if isempty(options.outerRadius),options.outerRadius=3*truth.totalLength;end
validateattributes(options.kernelX,{'numeric'},{'real','finite','scalar','>',0,'<',truth.totalLength});
validateattributes(options.outerRadius,{'numeric'},{'real','finite','scalar','positive'});
c=options.kernelX;rel=original-[c 0];originalTheta=atan2(rel(:,2),rel(:,1));
assert(originalTheta(1)==0&&originalTheta(end)==pi&&all(diff(originalTheta)>0&diff(originalTheta)<pi), ...
    'RocketDomain:Kernel','给定kernel不在整个真实鼻/收缩/后体/底面的严格星形核中。');
kernelDistance=(rel(1:end-1,1).*rel(2:end,2)-rel(1:end-1,2).*rel(2:end,1))./vecnorm(diff(original),2,2);
assert(all(isfinite(kernelDistance)&kernelDistance>0),'RocketDomain:Kernel','内核至每个真实墙段须为严格正距离。');
wall=original(1,:);parentSegment=[];wallRegions={};
for i=1:size(original,1)-1
    if strcmp(originalRegions{i},'base'),pins=linspace(0,1,options.baseRadialSegments+1);else,pins=[0 1];end
    A=original(i,:);E=original(i+1,:)-A;
    for k=1:numel(pins)-1
        left=A+pins(k)*E;right=A+pins(k+1)*E;
        ta=atan2(left(2),left(1)-c);tb=atan2(right(2),right(1)-c);
        n=max(1,ceil(rad2deg(tb-ta)/options.maxWallAngleStepDeg));angles=linspace(ta,tb,n+1);
        for h=2:n+1
            if h==n+1,point=right;
            else
                e=[cos(angles(h)),sin(angles(h))];d=A-[c 0];den=e(1)*E(2)-e(2)*E(1);
                t=(d(1)*e(2)-d(2)*e(1))/den;
                assert(isfinite(t)&&t>pins(k)&&t<pins(k+1),'RocketDomain:RayIntersection','射线与原墙段交点无效。');
                point=A+t*E;
            end
            wall(end+1,:)=point;parentSegment(end+1,1)=i;wallRegions(end+1,1)=originalRegions(i); %#ok<AGROW>
        end
    end
end
theta=atan2(wall(:,2),wall(:,1)-c);rhoWall=hypot(wall(:,1)-c,wall(:,2));
assert(all(diff(theta)>0&diff(theta)<pi)&&all(rhoWall>0),'RocketDomain:Kernel','加点后射线顺序错误。');
span=options.outerRadius-rhoWall;nLayers=options.radialLayers;first=options.firstRadialSpacing;
assert(all(isfinite(span)&span>0)&&first*nLayers<=min(span), ...
    'RocketDomain:OuterBoundary','外域不足或首径向间距过大，无法建立向外非减增长的正层。');
X=zeros(size(wall,1),nLayers+1);R=X;growth=zeros(size(wall,1),1);
for i=1:size(wall,1)
    growth(i)=growthRatio(first,span(i),nLayers);
    increments=first*growth(i).^(0:nLayers-1);rho=rhoWall(i)+[0 cumsum(increments)];rho(end)=options.outerRadius;
    assert(all(isfinite(rho))&&all(diff(rho)>0),'RocketDomain:RadialLayers','径向层必须严格递增。');
    X(i,:)=c+rho*cos(theta(i));R(i,:)=rho*sin(theta(i));
end
% 真正轴边和原墙端点直接保留，避免sin(pi)残差或重算几何造成伪间隙。
X(1,:)=c+hypot(X(1,:)-c,R(1,:));R(1,:)=0;
X(end,:)=c-hypot(X(end,:)-c,R(end,:));R(end,:)=0;
X(:,1)=wall(:,1);R(:,1)=wall(:,2);
fv=buildAxisymmetricFVGeometry(X,R);

%% 三、逐墙段绑定与完整区域账；映射候选不冒充已完成场量投影
nWall=size(wall,1)-1;wallIDs=find(fv.faces.boundarySide==3);order=zeros(nWall,1);
for i=1:nWall
    owner=sub2ind(fv.shape,i,1);id=wallIDs(fv.faces.owner(wallIDs)==owner);
    assert(isscalar(id),'RocketDomain:WallBinding','每个真实墙段须恰好对应一个边界面。');order(i)=id;
end
wallIDs=order;bands=originalBands(parentSegment);candidates=cell(nWall,1);
for i=1:nWall,candidates{i}=find(p.faceBandID(:)==bands(i));end
lengths=vecnorm(diff(wall),2,2);wallArea=.5*(wall(1:end-1,2)+wall(2:end,2)).*lengths;
originalArea=.5*(original(1:end-1,2)+original(2:end,2)).*vecnorm(diff(original),2,2);
grouped=accumarray(parentSegment,wallArea,[numel(originalArea) 1]);
assert(max(abs(grouped-originalArea)./max(originalArea,realmin))<2e-12&& ...
    all(abs(fv.faces.areaPerRadian(wallIDs)-wallArea)<=512*eps*max(wallArea,realmin)), ...
    'RocketDomain:WallArea','墙段加点不能漏面积、归一化q或改变原母线。');
axisIDs={find(fv.faces.boundarySide==1),find(fv.faces.boundarySide==2)};
allAxis=[axisIDs{1};axisIDs{2}];outerIDs=find(fv.faces.boundarySide==4);
assert(all(fv.faces.areaPerRadian(allAxis)==0)&&all(fv.faces.centroid(allAxis,2)==0)&& ...
    all(fv.faces.areaPerRadian(outerIDs)>0),'RocketDomain:Axis','轴必须是真r=0且积分面积严格为0。');
bodyVolume=frustumVolumePerRadian(original);outer=[X(:,end),R(:,end)];outerVolume=frustumVolumePerRadian(outer);
volume=sum(fv.volumePerRadian);expected=outerVolume-bodyVolume;
assert(abs(volume-expected)<=2e-12*(abs(volume)+abs(expected)), ...
    'RocketDomain:Volume','全域r权体积与真实外边界减实体体积不符。');
base=strcmp(wallRegions,'base');baseArea=2*pi*sum(wallArea(base));
assert(abs(baseArea-pi*truth.baseRadius^2)<=1024*eps*pi*truth.baseRadius^2, ...
    'RocketDomain:Base','真实底面必须从基底外缘完整到轴中心，不得加虚盖或漏环。');
axisSkew=fv.faces.centroid(allAxis,1)-fv.volumeCentroid(fv.faces.owner(allAxis),1);
domain=struct('kind','complete-native-rocket-star-ray-exterior-domain-v1','options',options, ...
    'coordinateConvention','SI; nose +X; zero-angle freestream -X; full per-radian FV', ...
    'wallPolylineXR',wall,'wallEndpointsXR',cat(3,wall(1:end-1,:),wall(2:end,:)), ...
    'wallBoundaryFaceID',wallIDs,'wallOwnerCell',fv.faces.owner(wallIDs),'wallSourceFaceBandID',bands, ...
    'wallSourceOriginalSegment',parentSegment,'wallRegion',{wallRegions},'wallSourceTriangleIDs',{candidates}, ...
    'sourceTriangleMappingKind','all angular triangles of matching original band are projection candidates; no field projection performed', ...
    'originalWallPolylineXR',original,'originalWallVertexID',originalVertexIDs,'originalWallFaceBandID',originalBands, ...
    'originalWallRegion',{originalRegions},'sourceVertices',V,'sourceFaces',F,'sourceFaceBandID',p.faceBandID(:), ...
    'kernel',[c 0],'kernelDistanceToOriginalWallM',kernelDistance,'wallPolarAngleRad',theta, ...
    'radialGrowthRatio',growth,'axisBoundarySide',[1 2],'axisBoundaryFaceID',{axisIDs}, ...
    'axisFaceCentroidXR',fv.faces.centroid(allAxis,:),'axisFaceOwnerXOffsetM',axisSkew, ...
    'requiredAxisGradientPolicy','mirror-owner-parity-v1; legacy AxisSkew must not be bypassed', ...
    'outerBoundaryFaceID',outerIDs,'outerBoundaryPolylineXR',outer, ...
    'wallAreaPerRadianM2',wallArea,'fullRevolvedWallAreaM2',2*pi*sum(wallArea), ...
    'sourceTriangulatedWallAreaM2',sum(sourceAreas),'analyticBodyWallAreaM2',truth.exactBaselineArea, ...
    'fullBaseAreaM2',baseArea,'bodyVolumePerRadianM3',bodyVolume,'sourceTriangulatedBodyVolumeM3',sourceVolume, ...
    'analyticBodyVolumeM3',truth.exactBaselineVolume,'outerPolygonVolumePerRadianM3',outerVolume, ...
    'fluidVolumePerRadianM3',volume,'wallNormalDistanceM',fv.faces.normalDistance(wallIDs), ...
    'firstRadialSpacingIsWallNormalSpacing',false,'nearWallResolutionValidated',false,'farfieldAdequacyValidated',false, ...
    'metricClosureRelative',fv.metricClosureRelative,'endcapsAdded',false,'physicalSolution',false, ...
    'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
end

function q=growthRatio(first,span,count)
% 模块用途：仅求几何层厚增长率；端点强制命中给定远场，不裁任何流动状态。
if abs(span-first*count)<=128*eps*max(span,first*count),q=1;return;end
lo=1;hi=2;
while first*sum(hi.^(0:count-1))<span
    hi=2*hi;assert(isfinite(hi),'RocketDomain:RadialLayers','径向增长率不可表示。');
end
for k=1:80
    mid=(lo+hi)/2;
    if first*sum(mid.^(0:count-1))<span,lo=mid;else,hi=mid;end
end
q=(lo+hi)/2;
end
function volume=frustumVolumePerRadian(points)
% 模块用途：绕X轴旋转真实折线的独立圆台体积；竖直真实底面不产生额外体积。
assert(all(diff(points(:,1))<=512*eps*max(1,max(abs(points(:,1))))), ...
    'RocketDomain:Volume','当前完整箭体/外圆边界须沿X单调不增。');
r=points(:,2);volume=sum(-diff(points(:,1)).*(r(1:end-1).^2+r(1:end-1).*r(2:end)+r(2:end).^2))/6;
end
