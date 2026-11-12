function estimate = estimateNoseRadius(mesh)
%ESTIMATENOSERADIUS 从已对齐的三维表面网格识别轴对称球形钝鼻半径。
% 模块作用：为 Sutton-Graves 驻点公式提供可审查的几何半径来源。
% 在鼻尖起算的三个尺度分别拟合球面，并检查拟合残差、半径一致性、
% 球心与 +X 轴的一致性，以及几何最前端是否与拟合球极点重合。
% 尖锥、平封盖、采样过稀或非球形鼻部返回 available=false。
% 该半径是 STL 几何半径，不验证激波、壁面化学或 Sutton-Graves 适用性。

%% 1. 按 X 截面确定鼻端点及三个拟合尺度
V=double(mesh.V);
validateattributes(V,{'numeric'},{'real','finite','ncols',3});
estimate=struct('available',false,'radius',NaN, ...
    'method','鼻端多尺度球面拟合与轴线/极点一致性检查', ...
    'scope','只识别轴对称球形钝鼻的 STL 几何半径；不代表流动求解结果', ...
    'status','几何条件尚未检查', ...
    'fitRadii',nan(1,3),'fitRelativeRMS',nan(1,3), ...
    'axisOffsetFraction',NaN,'tipMismatchFraction',NaN, ...
    'radiusSpreadFraction',NaN,'outerDepthRadiusFraction',NaN, ...
    'fitSectionCount',0,'minimumSectionCount',6);
if size(V,1)<50
    estimate.status='鼻端节点过少'; return;
end
xTip=max(V(:,1)); L=xTip-min(V(:,1));
if L<=0, estimate.status='X 向长度为零'; return; end
tol=max(L*1e-8,1e-10);
depth=xTip-V(:,1);
tipIDs=find(depth<=tol);
if numel(tipIDs)>3
    estimate.status='最前端为多节点平端或无法确定单一球极点'; return;
end
tip=mean(V(tipIDs,:),1);
levels=unique(round(depth/tol)*tol);
levels=levels(levels>tol);
if numel(levels)<6
    estimate.status='鼻端轴向截面不足，无法验证多尺度半径'; return;
end
%% 2. 各尺度独立球面拟合，记录几何质量和参数稳定性
% 从较大鼻端前缀向内检查；固定取 24 圈会把半球后的柱段混入拟合。
% 每个候选仍需三个嵌套尺度、至少六圈及下方全部几何判据。
% 不接受孤立三点球拟合；优先保留最宽的已验证球形范围。
for outer=min(24,numel(levels)):-1:6
levelIndices=unique(round([1/3 2/3 1]*outer));
fits=cell(1,3);
for j=1:3
    selected=depth<=levels(levelIndices(j))+tol;
    fits{j}=fitSphere(V(selected,:));
    estimate.fitRadii(j)=fits{j}.radius;
    estimate.fitRelativeRMS(j)=fits{j}.relativeRMS;
    if ~fits{j}.valid
        break;
    end
end
if any(cellfun(@isempty,fits)) || ~all(cellfun(@(f)f.valid,fits))
    continue;
end
R=median(estimate.fitRadii);
outerFit=fits{3};
estimate.radiusSpreadFraction= ...
    (max(estimate.fitRadii)-min(estimate.fitRadii))/R;
estimate.axisOffsetFraction=norm(outerFit.center(2:3)-tip(2:3))/R;
estimate.tipMismatchFraction=abs(outerFit.center(1)+outerFit.radius-xTip)/R;
estimate.outerDepthRadiusFraction=levels(levelIndices(3))/R;
estimate.fitSectionCount=outer;

%% 3. 保守地批准球形钝鼻；未通过时不输出用于热流的半径
valid=all(estimate.fitRelativeRMS<=0.003) && ...
    estimate.radiusSpreadFraction<=0.08 && ...
    estimate.axisOffsetFraction<=0.03 && ...
    estimate.tipMismatchFraction<=0.02 && ...
    estimate.outerDepthRadiusFraction>=0.08 && ...
    estimate.outerDepthRadiusFraction<=2.05;
if valid
    estimate.available=true;
    estimate.radius=R;
    estimate.status='多尺度球形钝鼻几何检查通过';
    return;
end
end
estimate.status='球形钝鼻几何一致性不足；驻点半径未自动采用';
end

function fit=fitSphere(points)
% 局部球面最小二乘：中心化缩放后解线性代数球，再计算真实径向残差。
fit=struct('valid',false,'radius',NaN,'center',[NaN NaN NaN], ...
    'relativeRMS',NaN);
if size(points,1)<20, return; end
center=mean(points,1);
scale=max(vecnorm(points-center,2,2));
if scale<=0, return; end
q=(points-center)/scale;
A=[2*q,ones(size(q,1),1)];
singular=svd(A,0);
if singular(end)<1e-8*singular(1), return; end
z=A\sum(q.^2,2);
radiusSquared=sum(z(1:3).^2)+z(4);
if radiusSquared<=0, return; end
fit.center=center+scale*z(1:3)';
fit.radius=scale*sqrt(radiusSquared);
residual=vecnorm(points-fit.center,2,2)-fit.radius;
fit.relativeRMS=sqrt(mean(residual.^2))/fit.radius;
fit.valid=isfinite(fit.relativeRMS) && fit.radius>0;
end
