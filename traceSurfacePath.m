function result = traceSurfacePath(mesh,faceID)
%TRACESURFACEPATH 重建指定面元的三维离散准流线路径。
% 模块作用：从 traceSurfaceStreamlines 的前驱树回溯到上游种子，
% 将种子面心、共享边中点和各下游面心连成可绘制的三维折线。
% 折线总长与供边界层积分的离散流动距离一致；不是连续 CFD 流线。

%% 1. 检查面元和已计算的路径树
nf=size(mesh.F,1);
validateattributes(faceID,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1,'<=',nf});
result=struct('available',false, ...
    'method','共享边中点与面元前驱树的三维折线重建', ...
    'scope','网格上的离散准流线；不表示连续流场或分离再附', ...
    'status','所选面元没有有效流线路径', ...
    'faceID',faceID,'faceIDs',zeros(0,1), ...
    'points',zeros(0,3),'distance',NaN);
if ~isfield(mesh,'surfaceStreamlines') || ...
        ~isfield(mesh.surfaceStreamlines,'crossingPoint') || ...
        ~isfield(mesh,'localGeometry')
    result.status='缺少含共享边穿越点的流线结果';
    return;
end
s=mesh.surfaceStreamlines;
if ~s.connectedFace(faceID)
    return;
end

%% 2. 沿前驱链回溯；显式检查环，避免坏路径无限循环
chain=zeros(nf,1);
visited=false(nf,1);
n=0; id=faceID;
while id>0
    if visited(id)
        result.status='路径前驱图包含环';
        return;
    end
    visited(id)=true;
    n=n+1; chain(n)=id;
    id=s.predecessorFace(id);
end
chain=flipud(chain(1:n));
C=mesh.localGeometry.center;
points=zeros(2*n-1,3);
points(1,:)=C(chain(1),:);
for k=2:n
    point=s.crossingPoint(chain(k),:);
    if any(~isfinite(point))
        result.status='前驱链缺少共享边穿越点';
        return;
    end
    points(2*k-2,:)=point;
    points(2*k-1,:)=C(chain(k),:);
end
distance=sum(vecnorm(diff(points,1,1),2,2));
if abs(distance-s.distance(faceID))>1e-7*max(1,distance)
    result.status='重建折线与流动距离不一致';
    return;
end

%% 3. 保存路径面元、空间点和真实折线长度
result.available=true;
result.status='三维离散准流线已重建';
result.faceIDs=chain;
result.points=points;
result.distance=distance;
end
