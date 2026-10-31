function result = propagateSurfaceSeparation(path,onsetFace)
%PROPAGATESURFACESEPARATION 沿离散准流线传递疑似层流分离的无解状态。
% 模块作用：已有 Thwaites 判据发现层流分离后，沿 Dijkstra 前驱树
% 标记全部下游面元。当前程序没有分离泡或再附模型，不能在下一面元
% 因局部速度梯度恢复而自动重新计算附着流热流。
% 传递仅是适用性筛选，不计算分离区热流，也不预测真实再附位置。

%% 1. 校验路径与分离起点，防止旧网格或循环路径产生假结果
required={'predecessorFace','pathCost','connectedFace'};
for k=1:numel(required)
    if ~isfield(path,required{k})
        error('RocketAeroHeat:Separation:Path','路径缺少 %s。',required{k});
    end
end
parent=double(path.predecessorFace(:));
cost=double(path.pathCost(:));
connected=logical(path.connectedFace(:));
onset=logical(onsetFace(:));
nf=numel(parent);
if numel(cost)~=nf || numel(connected)~=nf || numel(onset)~=nf || ...
        any(parent<0 | parent>nf | parent~=floor(parent)) || ...
        any(onset & ~connected) || any(~isfinite(cost(connected)))
    error('RocketAeroHeat:Separation:Path','分离标志与有效前驱路径不一致。');
end

%% 2. 按路径成本递增传递；一旦分离，不自动假设再附
blocked=onset;
ids=find(connected);
[~,order]=sort(cost(ids));
ids=ids(order);
for k=1:numel(ids)
    id=ids(k); upstream=parent(id);
    if upstream==0, continue; end
    if ~connected(upstream) || cost(upstream)>=cost(id)
        error('RocketAeroHeat:Separation:Path', ...
            '前驱必须是成本更小的有效上游面元。');
    end
    blocked(id)=blocked(id) || blocked(upstream);
end

%% 3. 保存来源与可用于热流空白诊断的面元标志
result=struct;
result.method='Thwaites 层流分离判据 + 三维准流线前驱树下游传递';
result.scope=['仅筛除缺少分离/再附模型的路径；' ...
    'lambda<=-0.09 为疑似层流分离，不是三维分离线的精确位置'];
result.onsetFace=onset;
result.unresolvedFace=blocked;
result.downstreamFace=blocked & ~onset;
result.summary=struct('onsetFaces',nnz(onset), ...
    'unresolvedFaces',nnz(blocked), ...
    'downstreamFaces',nnz(blocked & ~onset));
end
