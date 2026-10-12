function result = calcBluntNoseHeating(mesh,freestream,wallTemperature)
%CALCBLUNTNOSEHEATING 将驻点单点热流延拓到已验证的球形鼻端面元。
% 模块作用：只在 STL 多尺度球拟合通过且驻点热流可用时，
% 用相对逆来流方向的球面极角余弦权重填补鼻端面元。
% 换热系数同样按 cos(theta) 缩放，恢复温度采用驻点总温。
% 这是一种以 Sutton-Graves 驻点值为锚点的低阶空间延拓，
% 不是 Fay-Riddell 边界层解，也未验证激波层沿球面热流分布。
% 球鼻允许合成入射角 <=30°，仍只输出几何前半球中相对偏置驻点
% theta<=80° 的面元；其他面元保持 NaN。
% 可选逐面元壁温使球鼻焓差逐面元更新；驻点单点仍保留自己的
% 壁温和参考热流，不用某一鼻部面元温度替代整个鼻端温度。

%% 1. 初始化与严格检查适用范围
nf=size(mesh.F,1);
if nargin>=3
    validateattributes(wallTemperature,{'numeric'}, ...
        {'real','vector','numel',nf});
    wallTemperature=double(wallTemperature(:));
    if any(isinf(wallTemperature) | wallTemperature<=0)
        error('RocketAeroHeat:NoseHeating:WallTemperature','面元壁温须为正数或 NaN。');
    end
end
result=struct('available',false, ...
    'method','Sutton-Graves 偏置驻点值 + 球面迎流极角余弦几何延拓', ...
    'scope',['已验证球形钝鼻、合成入射角 <=30°、前半球且' ...
        '相对偏置驻点 theta<=80°；未经局部激波层/边界层验证'], ...
    'status','缺少球形鼻端或有效驻点热流', ...
    'heatFlux',nan(nf,1),'heatTransferCoefficient',nan(nf,1), ...
    'recoveryTemperature',nan(nf,1), ...
    'faceMask',false(nf,1),'polarAngleDeg',nan(nf,1), ...
    'axisPolarAngleDeg',nan(nf,1), ...
    'stagnationHeatFlux',NaN,'noseRadius',NaN, ...
    'stagnationPosition',[NaN NaN NaN]);
if ~isfield(mesh,'stagnation') || ...
        ~isfield(mesh.stagnation,'noseRadiusEstimate') || ...
        ~mesh.stagnation.noseRadiusEstimate.available || ...
        ~isfield(mesh,'stagnationHeating') || ...
        ~mesh.stagnationHeating.available || ...
        ~isfield(mesh,'localGeometry')
    return;
end
if hypot(freestream.alpha,freestream.beta)>30
    result.status='合成入射角超过 30°，偏置球鼻延拓不适用';
    return;
end
R=mesh.stagnation.noseRadiusEstimate.radius;
if abs(mesh.stagnationHeating.radius/R-1)>0.05
    result.status='驻点热流使用的半径与 STL 球形鼻端不一致';
    return;
end

%% 2. 从真实 STL 节点检查面元确实位于球形前半部
V=double(mesh.V); F=double(mesh.F); g=mesh.localGeometry;
xTip=max(V(:,1));
tipIDs=abs(V(:,1)-xTip)<=max(1e-9,1e-8*R);
tipYZ=mean(V(tipIDs,2:3),1);
sphereCenter=[xTip-R,tipYZ];
radialVertices=vecnorm(V-sphereCenter,2,2);
vertexError=abs(radialVertices-R)/R;
faceError=max(vertexError(F),[],2);
radialCenters=g.center-sphereCenter;
centerLength=vecnorm(radialCenters,2,2);
cosAxis=radialCenters(:,1)./centerLength;
upstream=-double(freestream.unitFlowDirection(:)');
cosTheta=(radialCenters*upstream')./centerLength;
normalAlignment=sum(g.normal.*radialCenters,2)./centerLength;
front=cosAxis>=cosd(80) & cosAxis<=1+1e-6 & ...
    cosTheta>=cosd(80) & cosTheta<=1+1e-6;
vertexX=reshape(V(F(:),1),size(F));
mask=front & faceError<=0.03 & ...
    normalAlignment>=cosd(20) & ~g.topologyIssue & ...
    all(vertexX>=sphereCenter(1)-0.01*R,2);
if ~any(mask)
    result.status='没有通过球面位置与法向检查的前半球面元';
    return;
end

%% 3. 驻点锚定的角度延拓；未覆盖面元保持 NaN
weight=max(0,min(1,cosTheta(mask)));
stag=mesh.stagnationHeating;
result.heatFlux(mask)=stag.heatFlux*weight;
result.heatTransferCoefficient(mask)= ...
    stag.heatTransferCoefficient*weight;
result.recoveryTemperature(mask)=stag.adiabaticWallTemperature;
if nargin>=3
    % 先由原驻点 h 消去原壁温的等效比热，再乘当前面元的等效比热。
    % Taw=Tw 时仍使用有限的 cp 极限，不通过 q/(Taw-Tw) 构造 0/0。
    factor=ones(nnz(mask),1);
    if isfield(stag,'thermalPropertyMode') && strcmp(stag.thermalPropertyMode,'thermallyPerfect')
        factor=airEnthalpySecantCp(stag.adiabaticWallTemperature,wallTemperature(mask)) / ...
            airEnthalpySecantCp(stag.adiabaticWallTemperature,stag.wallTemperature);
    end
    result.heatTransferCoefficient(mask)=result.heatTransferCoefficient(mask).*factor;
    result.heatFlux(mask)=result.heatTransferCoefficient(mask).* ...
        (stag.adiabaticWallTemperature-wallTemperature(mask));
end
result.polarAngleDeg(mask)=acosd(weight);
result.axisPolarAngleDeg(mask)=acosd(max(-1,min(1,cosAxis(mask))));
result.faceMask=mask;
% 几何覆盖不等于热流有解；未知/物性越界壁温不能留下有效标志。
result.faceMask=mask & isfinite(result.heatFlux) & ...
    isfinite(result.heatTransferCoefficient);
result.heatTransferCoefficient(~result.faceMask)=NaN;
result.recoveryTemperature(~result.faceMask)=NaN;
result.stagnationHeatFlux=stag.heatFlux;
result.noseRadius=R;
result.stagnationPosition=stag.position;
result.available=any(result.faceMask);
result.status=sprintf('球形鼻端低阶延拓覆盖 %d 个面元',nnz(result.faceMask));
end
