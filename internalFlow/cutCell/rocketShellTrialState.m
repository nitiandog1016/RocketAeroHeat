function candidate=rocketShellTrialState(Q,directionU,scale,alpha,gas,coordinates)
%ROCKETSHELLTRIALSTATE 为同一完整NS Newton切向方向构造显式非线性候选路径。
% 模块一：conservative保持原Q+alpha*dQ逐值表达式，默认数值路径不改。
% 模块二：primitive-tangent用链式法则得到drho/du/dT，更新后恢复全部五守恒量。
% 模块三：本函数不裁密度/温度/内能，不推进物理时间；候选仍须原完整残差验收。
assert(ismember(coordinates,{'conservative','primitive-tangent'}),'ShellTrial:Coordinates','未知候选状态坐标。');
validateattributes(alpha,{'numeric'},{'scalar','real','finite','nonnegative'});
scale=reshape(scale,1,5);validateattributes(scale,{'numeric'},{'numel',5,'finite','real','positive'});
assert(size(Q,2)==5&&isequal(size(Q),size(directionU))&&isreal(Q)&&isreal(directionU)&& ...
    all(isfinite(Q),'all')&&all(isfinite(directionU),'all'),'ShellTrial:Shape','候选切向方向/完整五量不符。');
if strcmp(coordinates,'conservative')||alpha==0
    candidate=Q+alpha*directionU.*scale;return;
end
assert(isstruct(gas)&&isscalar(gas.R)&&isfinite(gas.R)&&gas.R>0&&isscalar(gas.gamma)&& ...
    isfinite(gas.gamma)&&gas.gamma>1,'ShellTrial:Gas','需要明确的定比热气体。');
rho=Q(:,1);assert(all(rho>0),'CutNS:State','基态密度非正，不裁剪。');
cv=gas.R/(gas.gamma-1);u=Q(:,2:4)./rho;
T=(Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);assert(all(T>0),'CutNS:InternalEnergy','基态内能非正。');
dQ=directionU.*scale;drho=dQ(:,1);du=(dQ(:,2:4)-u.*drho)./rho;
dT=(dQ(:,5)-sum(u.*dQ(:,2:4),2)+(.5*sum(u.^2,2)-cv*T).*drho)./(rho*cv);
newRho=rho+alpha*drho;newU=u+alpha*du;newT=T+alpha*dT;
candidate=[newRho,newRho.*newU,newRho.*(cv*newT+.5*sum(newU.^2,2))];
end
