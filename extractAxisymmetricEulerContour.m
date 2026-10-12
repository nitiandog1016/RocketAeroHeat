function [contour,geometry] = extractAxisymmetricEulerContour(mesh)
%EXTRACTAXISYMMETRICEULERCONTOUR 从原生截面拓扑提取回转体母线。
% 模块作用：检验圆度、中心线与尖鼻条件，拒绝直接把有翼/偏心模型
% 平均成轴对称体。容许 2% 半径的 STL/重采样离散误差。
if ~isfield(mesh,'param') || ~isfield(mesh.param,'sectionNodeIDs') || ...
        ~isfield(mesh.param,'sectionX')
    error('RocketAeroHeat:Euler:Topology','需要 autoRemeshSurface 的原生截面拓扑。');
end
x=mesh.param.sectionX(:); ids=mesh.param.sectionNodeIDs;
n=numel(x); centers=nan(n,2); radii=zeros(n,1); errorRadius=zeros(n,1);
for j=1:n
    yz=mesh.V(ids{j},2:3);
    if size(yz,1)==1, centers(j,:)=yz; continue;end
    if size(yz,1)<3
        error('RocketAeroHeat:Euler:Roundness','截面节点不足以检验圆度。');
    end
    % 平移后拟合圆心，避免较大模型坐标导致法方程精度损失。
    base=mean(yz,1); q=yz-base; A=[2*q,ones(size(q,1),1)];
    if rank(A)<3,error('RocketAeroHeat:Euler:Roundness','截面退化。');end
    fit=A\sum(q.^2,2); centers(j,:)=base+fit(1:2).';
    radius=vecnorm(yz-centers(j,:),2,2); radii(j)=mean(radius);
    errorRadius(j)=max(abs(radius-radii(j)));
end
[x,order]=sort(x); centers=centers(order,:); radii=radii(order); errorRadius=errorRadius(order);
% 重构器把鼻尖作为封口节点存储，不属于普通截面环。
if isfield(mesh.param,'noseCapNodeID') && ~isempty(mesh.param.noseCapNodeID) && ...
        isfield(mesh.param,'noseCapType') && strcmp(mesh.param.noseCapType,'tip')
    tip=mesh.V(mesh.param.noseCapNodeID,:);
    if size(tip,1)==1 && tip(1)>x(end)
        x(end+1)=tip(1);centers(end+1,:)=tip(2:3);radii(end+1)=0;errorRadius(end+1)=0;
    end
end
maxRadius=max(radii); center=median(centers,1);
if maxRadius<=0 || any(diff(x)<=0) || radii(end)>1e-6*maxRadius || ...
        any(radii(1:end-1)<=0)
    error('RocketAeroHeat:Euler:Shape','要求单一尖鼻、正半径侧壁和不同的截面 X。');
end
centerError=max(vecnorm(centers-center,2,2))/maxRadius;
roundError=max(errorRadius./max(radii,0.01*maxRadius));
if centerError>0.02 || roundError>0.02
    error('RocketAeroHeat:Euler:Axisymmetry', ...
        '轴对称检验未通过：中心线偏差 %.3g，圆度误差 %.3g；阈值均为 0.02。',centerError,roundError);
end
radii(end)=0; contour=[x,radii];
geometry=struct('centerYZ',center,'relativeCenterError',centerError, ...
    'relativeRoundnessError',roundError,'method','截面圆拟合 + 中心线一致性检验', ...
    'scope','允许 2% 离散几何误差；母线分片线性，不含尾盖/尾迹');
end
