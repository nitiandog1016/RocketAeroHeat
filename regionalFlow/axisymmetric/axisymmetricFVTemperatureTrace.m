function trace=axisymmetricFVTemperatureTrace(mesh,T,gx,gr,boundaryT)
%AXISYMMETRICFVTEMPERATURETRACE 对本构系数采样温度作逐单元Barth-Jespersen重构。
% 模块一：温度上下界只来自本单元、面邻居和实际边界迹，不用全局温域裁温度。
% 模块二：每单元统一alpha缩放到所有面及平面源点的温度增量，保留单元T/Q。
% 模块三：只服务mu/k评估；原WLS热/黏性梯度及守恒通量结构不改变。
% 这不是完整NS状态正性限制器；alpha/原始迹/局部界显式保存用于独立复算。
nc=mesh.cellCount;f=mesh.faces;owner=f.owner;other=f.neighbour;inside=other>0;bf=mesh.boundaryFaces;
assert(isequal(size(T),[nc 1])&&isequal(size(gx),size(T))&&isequal(size(gr),size(T))&& ...
    isequal(size(boundaryT),[numel(bf) 1])&&all(isfinite([T;gx;gr;boundaryT]))&&all([T;boundaryT]>0), ...
    'AxisFV:TemperatureTrace','重构需要完整有限正单元/边界温度和真实梯度。');
% 按每条真实邻接批量归并，避免每个Jacobian扰动再逐面动态更新上下界。
boundCell=[(1:nc)';owner(inside);other(inside);owner(bf)];
boundTemperature=[T;T(other(inside));T(owner(inside));boundaryT];
lower=accumarray(boundCell,boundTemperature,[nc 1],@min);
upper=accumarray(boundCell,boundTemperature,[nc 1],@max);
leftDisplacement=f.centroid-mesh.volumeCentroid(owner,:);
leftDelta=gx(owner).*leftDisplacement(:,1)+gr(owner).*leftDisplacement(:,2);
rightDisplacement=f.centroid(inside,:)-f.neighbourShift(inside,:)-mesh.volumeCentroid(other(inside),:);
rightDelta=gx(other(inside)).*rightDisplacement(:,1)+gr(other(inside)).*rightDisplacement(:,2);
sourceDisplacement=mesh.planarCentroid-mesh.volumeCentroid;
sourceDelta=gx.*sourceDisplacement(:,1)+gr.*sourceDisplacement(:,2);
targetCell=[owner;other(inside);(1:nc)'];increments=[leftDelta;rightDelta;sourceDelta];
candidates=ones(size(increments));positive=increments>0;negative=increments<0;
candidates(positive)=(upper(targetCell(positive))-T(targetCell(positive)))./increments(positive);
candidates(negative)=(lower(targetCell(negative))-T(targetCell(negative)))./increments(negative);
% 极小非零斜率可给出+Inf比例；min(1,ratio)=1，不把有限输入误拒绝。
assert(all(~isnan(candidates)&candidates>=0),'AxisFV:TemperatureTrace','局部单调重构比例无效。');
alpha=accumarray(targetCell,min(1,candidates),[nc 1],@min,1);
% 只在已限制的斜率上留浮点内侧余量，避免边界值因一次加法向外舍入。
limited=alpha<1;alpha(limited)=alpha(limited)*(1-16*eps);
left=T(owner)+alpha(owner).*leftDelta;
right=left;right(inside)=T(other(inside))+alpha(other(inside)).*rightDelta;right(bf)=boundaryT;
faceTemperature=left;faceTemperature(inside)=.5*(left(inside)+right(inside));faceTemperature(bf)=boundaryT;
sourceTemperature=T+alpha.*sourceDelta;
rawFace=T(owner)+leftDelta;
rawFace(inside)=.5*(T(owner(inside))+leftDelta(inside)+T(other(inside))+rightDelta);
rawFace(bf)=boundaryT;
rawLeft=T(owner)+leftDelta;rawRight=rawLeft;
rawRight(inside)=T(other(inside))+rightDelta;rawRight(bf)=boundaryT;
trace=struct('kind','cellwise-Barth-Jespersen-constitutive-temperature-v1', ...
    'cellLower',lower,'cellUpper',upper,'cellLimiter',alpha,'limitedCellCount',nnz(limited), ...
    'faceTemperature',faceTemperature,'sourceTemperature',sourceTemperature, ...
    'rawFaceTemperature',rawFace,'rawSourceTemperature',T+sourceDelta, ...
    'rawLeftTemperature',rawLeft,'rawRightTemperature',rawRight, ...
    'leftTemperature',left,'rightTemperature',right,'cellStateClipped',false, ...
    'originalViscousThermalGradientLimited',false, ...
    'scope','Only constitutive mu/k sampling; not complete NS positivity or turbulence validation');
end
