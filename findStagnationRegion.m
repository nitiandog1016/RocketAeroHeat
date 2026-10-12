function mesh = findStagnationRegion(mesh,freestream,cfg)
%FINDSTAGNATIONREGION 寻找自由来流下的几何驻点候选区，供气动热准备。
% 调用：mesh=findStagnationRegion(mesh,freestream,cfg);
% 输入：最终网格、localGeometry、当前来流；不要求压力已计算。
% 输出：mesh.stagnation。所有“驻点”均为候选，不是外流求解结果。
%
% 方法：沿来流反方向寻找上游区域，再用外法向迎流程度筛选；
% 补充最上游顶点邻接面，避免尖鼻法向不迎流时漏掉几何尖端。
% 对候选面做边连通分组，保留多个候选区，不将多驻点强并成一个。
% 不计算遮挡、激波、分离或真实表面流线；凹形体需要额外验证。
%
% cfg.stagnation 可选参数：
%   upstreamFraction=0.05：沿来流投影长度最前面的 5%%。
%   alignmentAngleDeg=20：正常迎流面与上游方向的夹角阈值。
%   peakMuTolerance=0.02：仅保留接近上游最大迎流程度的面。
%   physicalNoseRadius：用户提供的物理鼻部半径，m；没有则 NaN。
% 物理半径即使已提供，也不能单独证明某个热流模型适用。
%
% sphereFit 仅对候选区附近网格拟合球面作为诊断，其半径不会直接用于热流。
% 另由 estimateNoseRadius 在三个几何尺度验证球形钝鼻；验证通过时可
% 自动取得模型本身的几何鼻尖半径。尖锥和平封盖不会自动得到半径。

if nargin<3, cfg=struct; end
%% 1. 复用分类模块验证几何与来流的一致性，并使用未经阈值截断的 mu
mesh=classifyAerothermalRegion(mesh,freestream,cfg);
g=mesh.localGeometry; V=double(mesh.V); F=double(mesh.F);
nf=size(F,1); mu=mesh.aerothermalRegion.mu;
flow=mesh.aerothermalRegion.unitFlowDirection; upstream=-flow;
settings=struct('upstreamFraction',0.05,'alignmentAngleDeg',20,'peakMuTolerance',0.02);
if isfield(cfg,'stagnation')
    for key=fieldnames(settings)'
        if isfield(cfg.stagnation,key{1}), settings.(key{1})=cfg.stagnation.(key{1}); end
    end
end
validateattributes(settings.upstreamFraction,{'numeric'},{'real','finite','scalar','>',0,'<=',1});
validateattributes(settings.alignmentAngleDeg,{'numeric'},{'real','finite','scalar','>',0,'<',90});
validateattributes(settings.peakMuTolerance,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
% 使用相对坐标投影，避免模型平移后损失有效精度。
origin=mean(V,1); projection=(V-origin)*upstream';
front=max(projection); span=front-min(projection);
if span<=0, error('RocketAeroHeat:Stagnation:Extent','沿来流方向的几何长度必须非零。'); end
tol=max(span*1e-10,64*eps(max(abs(V(:)))));
depth=front-(g.center-origin)*upstream';
upstreamFace=depth<=settings.upstreamFraction*span;
leadingNodes=find(front-projection<=tol);
leadingFace=any(ismember(F,leadingNodes),2);
search=(upstreamFace|leadingFace) & mu>0;
if ~any(search)
    error('RocketAeroHeat:Stagnation:NoCandidate','上游区域未找到迎风面，请检查外法向。');
end

%% 2. 筛选迎流面并保留最前端邻接面，区分高迎流区域与纯几何前端
peakMu=max(mu(search)); normalLimit=cosd(settings.alignmentAngleDeg);
nearPeak=search & mu>=peakMu-settings.peakMuTolerance;
% 对整个迎流法向连续区域取连通分量，避免局部三角法向噪声把
% 同一钝鼻驻点区切碎为几十个孤立的“多驻点”。
strongAlignment=search & mu>=normalLimit;
candidate=strongAlignment | (leadingFace & mu>0);
if ~any(strongAlignment), candidate=candidate|nearPeak; end

%% 3. 共享边连通分组，同时构建局部拟合邻域
edges=sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2);
[~,~,group]=unique(edges,'rows'); count=accumarray(group,1);
ids=repmat((1:nf)',3,1);
first=accumarray(group,ids,[],@min); last=accumarray(group,ids,[],@max);
valid=count==2 & first~=last;
adj=sparse([first(valid);last(valid)],[last(valid);first(valid)],1,nf,nf);
candidateID=find(candidate);
components=conncomp(graph(adj(candidateID,candidateID)))';
componentID=zeros(nf,1); componentID(candidateID)=components;
componentArea=accumarray(components,g.area(candidateID));
% 数值法向造成的孤立单面元不应被解释成独立物理驻点。
significantComponent=componentArea>=0.05*sum(componentArea);
numberSignificant=nnz(significantComponent);
if numberSignificant==0
    significantComponent=true(size(componentArea));
    numberSignificant=numel(componentArea);
end
% 代表面优先取迎流程度最大者；并列时选择更靠上游者。
[~,order]=sortrows([-mu(candidateID),depth(candidateID)],[1 2]);
representative=candidateID(order(1));

%% 4. 在代表面附近扩展三圈邻域，做球面拟合诊断
patch=false(nf,1); patch(representative)=true;
for k=1:3, patch=patch|adj*double(patch)>0; end
nodes=unique(F(patch,:));
fit=fitSphereDiagnostic(V(nodes,:));
fit.faceIDs=find(patch); fit.nodeIDs=nodes;
fit.containsFeature=any(g.featureFace(patch));
fit.containsPoorAngle=any(g.poorAngle(patch));
% 球拟合即使残差很小，也不能自动证明热流半径有效。
fit.approvedForHeating=false;

%% 5. 整理候选信息和半径来源，不把推测值当成已知物理半径
s=struct;
s.version='1.0'; s.method='上游投影 + 法向迎流筛选 + 最前端邻接';
s.settings=settings; s.unitFlowDirection=flow;
s.upstreamFace=upstreamFace; s.leadingNodeIDs=leadingNodes;
s.leadingFace=leadingFace; s.candidateFace=candidate;
s.candidateFaceIDs=candidateID; s.componentID=componentID;
s.rawNumberComponents=max(components);
s.componentArea=componentArea;
s.significantComponent=significantComponent;
s.numberComponents=numberSignificant;
s.strongAlignmentFace=strongAlignment;
s.peakMu=peakMu; s.peakNormalAngleDeg=acosd(min(1,max(-1,peakMu)));
s.representativeFaceID=representative;
s.representativePosition=g.center(representative,:);
s.positionDefinition='候选代表面面心，不是精确驻点坐标';
s.candidateArea=sum(g.area(candidate));
s.sphereFit=fit;
s.noseRadiusEstimate=estimateNoseRadius(mesh);
s.physicalRadius=NaN; s.radiusSource='未提供；球形钝鼻几何检查未通过';
if s.noseRadiusEstimate.available
    s.physicalRadius=s.noseRadiusEstimate.radius;
    s.radiusSource='STL 鼻端多尺度球形钝鼻拟合';
end
if isfield(cfg,'stagnation') && isfield(cfg.stagnation,'physicalNoseRadius')
    radius=cfg.stagnation.physicalNoseRadius;
    validateattributes(radius,{'numeric'},{'real','finite','scalar','positive'});
    s.physicalRadius=radius; s.radiusSource='用户指定物理鼻部半径';
end
s.hasPhysicalRadius=isfinite(s.physicalRadius);
s.requiresFlowValidation=true; s.readyForHeating=false;
s.flags=struct('weakNormalAlignment',~any(strongAlignment), ...
    'multipleRegions',s.numberComponents>1,'featureCandidate',any(g.featureFace(candidate)), ...
    'poorAngleCandidate',any(g.poorAngle(candidate)), ...
    'topologyIssue',any(g.topologyIssue), ...
    'flatCapCandidate',any(g.isNoseCap(candidate)|g.isAftCap(candidate)));
if s.hasPhysicalRadius
    s.status='几何驻点候选与鼻尖半径已得到；仍须检查局部外流和经验式适用域';
else
    s.status='仅几何驻点候选；鼻尖半径未可靠确定';
end
mesh.stagnation=s;
end

function fit=fitSphereDiagnostic(points)
% 局部代数球拟合：中心化并缩放后求解，报告秩、条件数和径向残差。
% 平面点集秩不足；近似锥面也可能拟合出球，因此不能仅凭残差定形。
fit=struct('radius',NaN,'center',[NaN NaN NaN],'relativeRMS',NaN, ...
    'conditionNumber',Inf,'status','点数不足或局部退化');
if size(points,1)<10, return; end
center=mean(points,1); scale=max(vecnorm(points-center,2,2));
if scale<=0, return; end
q=(points-center)/scale; A=[2*q,ones(size(q,1),1)];
singular=svd(A,0); fit.conditionNumber=singular(1)/singular(end);
if singular(end)<1e-8*singular(1), fit.status='近似平面或秩不足'; return; end
solution=A\sum(q.^2,2); radiusSquared=sum(solution(1:3).^2)+solution(4);
if radiusSquared<=0, fit.status='球拟合半径无效'; return; end
fit.center=center+scale*solution(1:3)'; fit.radius=scale*sqrt(radiusSquared);
residual=vecnorm(points-fit.center,2,2)-fit.radius;
fit.relativeRMS=sqrt(mean(residual.^2))/fit.radius;
fit.status='未验证的局部球拟合；禁止直接用于热流公式';
end
