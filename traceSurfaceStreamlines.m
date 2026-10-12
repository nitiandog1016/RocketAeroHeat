function mesh = traceSurfaceStreamlines(mesh,cfg)
%TRACESURFACESTREAMLINES 根据三维外缘速度构建面元流向与流动距离。
% 算法：把速度投影到三角面切平面；共享边面元建立有向图；
% 以共享边中点近似跨面穿越点；两侧局部速度均须指向该穿越路径，
% 再以两段表面路径长度/流向一致性作权重；
% 从驻点候选区用 Dijkstra 求各面元的离散流动距离。
% 这仍是网格上的准流线近似，不是连续 CFD 流线积分；跨分离区不外推。
% 必须已有逐面元三维外缘速度及驻点候选区。

%% 1. 输入检查和速度切向投影
if nargin<2, cfg=struct; end
if ~isfield(mesh,'boundaryLayerEdge') || ...
        ~isfield(mesh.boundaryLayerEdge,'velocityVector')
    error('RocketAeroHeat:Streamlines:Velocity', ...
        '需要含 velocityVector 的逐面元三维外缘状态。');
end
if ~isfield(mesh,'stagnation') || ~isfield(mesh.stagnation,'candidateFace')
    error('RocketAeroHeat:Streamlines:Seeds','需要先确定驻点候选区。');
end
if ~isfield(mesh,'localGeometry')
    error('RocketAeroHeat:Streamlines:Geometry','缺少 localGeometry。');
end
F=double(mesh.F); V=double(mesh.V); nf=size(F,1);
g=mesh.localGeometry; C=double(g.center); N=double(g.normal);
U=double(mesh.boundaryLayerEdge.velocityVector);
validateattributes(C,{'numeric'},{'real','finite','size',[nf 3]});
validateattributes(N,{'numeric'},{'real','finite','size',[nf 3]});
validateattributes(U,{'numeric'},{'real','size',[nf 3]});
seeds=logical(mesh.stagnation.candidateFace(:));
if numel(seeds)~=nf || ~any(seeds)
    error('RocketAeroHeat:Streamlines:Seeds','驻点候选区为空或面元数不匹配。');
end
maxTurn=80;
if isfield(cfg,'streamline') && isfield(cfg.streamline,'maxStepDeviationDeg')
    maxTurn=cfg.streamline.maxStepDeviationDeg;
end
validateattributes(maxTurn,{'numeric'},{'real','finite','scalar','>',0,'<',90});
finiteVelocity=all(isfinite(U),2);
U(~finiteVelocity,:)=0;
Vt=U-sum(U.*N,2).*N;
speed=vecnorm(Vt,2,2);
valid=finiteVelocity & speed>max(1e-8,1e-8*max(speed));
direction=nan(nf,3);
direction(valid,:)=Vt(valid,:)./speed(valid);

%% 2. 共享边生成有向面元图，剔除逆流和大幅横向跳跃
edges=sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2);
[edgeVertices,~,group]=unique(edges,'rows');
counts=accumarray(group,1);
faceIDs=repmat((1:nf)',3,1);
first=accumarray(group,faceIDs,[],@min);
last=accumarray(group,faceIDs,[],@max);
shared=counts==2 & first~=last;
i=first(shared); j=last(shared);
mid=0.5*(V(edgeVertices(shared,1),:)+ ...
    V(edgeVertices(shared,2),:));
from=[i;j]; to=[j;i];
crossing=[mid;mid];
firstLeg=crossing-C(from,:);
secondLeg=C(to,:)-crossing;
firstLength=vecnorm(firstLeg,2,2);
secondLength=vecnorm(secondLeg,2,2);
surfaceStep=firstLength+secondLength;
alignmentFrom=sum(direction(from,:).*firstLeg,2)./firstLength;
alignmentTo=sum(direction(to,:).*secondLeg,2)./secondLength;
limit=cosd(maxTurn);
keep=valid(from) & valid(to) & isfinite(surfaceStep) & ...
    firstLength>0 & secondLength>0 & ...
    alignmentFrom>=limit & alignmentTo>=limit;
from=from(keep); to=to(keep);
weight=firstLength(keep)./alignmentFrom(keep)+ ...
    secondLength(keep)./alignmentTo(keep);

%% 3. 超级源点连接全部有效驻点候选面，单次 Dijkstra 求最短下游路径
source=nf+1;
seedIDs=find(seeds & valid);
seedSource='有效驻点候选面';
if isempty(seedIDs)
    % 驻点处真实切向速度为零：从与候选区相邻的第一圈有效面启动。
    boundary=(seeds(i) & valid(j));
    reverseBoundary=(seeds(j) & valid(i));
    seedIDs=unique([j(boundary);i(reverseBoundary)]);
    seedSource='驻点候选区相邻的第一圈有效面';
end
if isempty(seedIDs)
    % 脱体/尖鼻附近可能有数圈无附体解；使用最上游的有效外缘带。
    upstream=-mesh.stagnation.unitFlowDirection(:);
    projection=C*upstream;
    front=max(projection(valid));
    span=max(projection)-min(projection);
    seedIDs=find(valid & front-projection<=max(0.005*span,1e-9));
    seedSource='最上游的有效外缘带（驻点附近无附体解）';
end
if isempty(seedIDs)
    error('RocketAeroHeat:Streamlines:Seeds', ...
        '没有可用的外缘速度，无法追踪流线。');
end
G=digraph([from;repmat(source,numel(seedIDs),1)], ...
    [to;seedIDs],[weight;zeros(numel(seedIDs),1)],source);
pathCost=distances(G,source); pathCost=pathCost(1:nf)';
connected=isfinite(pathCost) & valid;
pathCost(~connected)=NaN;

%% 4. 沿最小成本路径累计真实几何长度，供边界层沿程积分
% Dijkstra 的权重 length/alignment 只用于选路，不能当作物理弧长。
% 对每个可达面元选取满足最短成本的前驱，再按成本递增累计步长。
predecessor=zeros(nf,1);
stepLength=nan(nf,1);
crossingPoint=nan(nf,3);
distance=nan(nf,1);
seedIDs=seedIDs(connected(seedIDs));
distance(seedIDs)=0;
stepLength(seedIDs)=0;
candidate=connected(from) & connected(to) & ~ismember(to,seedIDs);
cf=from(candidate); ct=to(candidate);
cl=surfaceStep(keep); cl=cl(candidate);
cp=crossing(keep,:); cp=cp(candidate,:);
cw=weight(candidate);
residual=abs(pathCost(cf)+cw-pathCost(ct));
% 排序后每个目标面元的首条边是最符合 Dijkstra 成本的前驱。
[~,order]=sortrows([ct,residual,cf],[1 2 3]);
ct=ct(order); cf=cf(order); cl=cl(order); cp=cp(order,:);
[uniqueTargets,firstIndex]=unique(ct,'stable');
predecessor(uniqueTargets)=cf(firstIndex);
stepLength(uniqueTargets)=cl(firstIndex);
crossingPoint(uniqueTargets,:)=cp(firstIndex,:);
[~,walkOrder]=sort(pathCost(connected));
connectedIDs=find(connected); connectedIDs=connectedIDs(walkOrder);
for k=1:numel(connectedIDs)
    id=connectedIDs(k);
    parent=predecessor(id);
    if parent>0 && isfinite(distance(parent))
        distance(id)=distance(parent)+stepLength(id);
    end
end
connected=connected & isfinite(distance);
pathCost(~connected)=NaN;
predecessor(~connected)=0;
stepLength(~connected)=NaN;
crossingPoint(~connected,:)=NaN;

%% 5. 保存方法、来源和覆盖度；不覆盖几何子午弧长
s=struct;
s.method='三维外缘速度投影 + 共享边两段有向图 Dijkstra 准流线';
s.scope='共享边中点近似穿越位置；未解析连续流线内插或分离再附';
s.velocitySource=mesh.boundaryLayerEdge.method;
s.distance=distance; s.pathCost=pathCost;
s.predecessorFace=predecessor; s.stepLength=stepLength;
s.crossingPoint=crossingPoint;
s.tangentDirection=direction;
s.tangentialSpeed=speed; s.connectedFace=connected;
s.seedFace=seeds; s.seedFaceIDs=seedIDs;
s.seedSource=seedSource;
s.maxStepDeviationDeg=maxTurn;
s.summary=struct('connectedFaces',nnz(connected), ...
    'unresolvedFaces',nnz(~connected),'seedFaces',numel(seedIDs));
mesh.surfaceStreamlines=s;
if isfield(mesh,'heating'), mesh=rmfield(mesh,'heating'); end
if isfield(mesh,'heatingCoverage'), mesh=rmfield(mesh,'heatingCoverage'); end
if isfield(mesh,'boundaryLayer'), mesh=rmfield(mesh,'boundaryLayer'); end
if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
end
