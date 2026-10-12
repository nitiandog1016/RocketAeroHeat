function result = calcWallConduction1D(h,Taw,time,material,initial,law)
%CALCWALLCONDUCTION1D 一维厚度导热，外侧对流/辐射与可选内侧换热。
% 默认背面绝热；material.backBoundary 启用时转入双面守恒导热核。
% 采用厚度方向有限体积离散、隐式 Euler 时间推进；
% 辐射 T^4 用 Newton 切线线性化，逐时刻求解三对角系统。
% 输入 h、Taw 为需一维求解的面元列向量；输出各面元表面时间历程。
% 可选 law 使 h 随表面温度迭代更新；外流及恢复温度不随时间变化。
% 本文件保留常物性快速分支；给定 propertyTable 时转入比热精确
% 积分/Kirchhoff 通量的变物性核；layers 输入转入逐层/接触热阻核。
% initial 可为统一初温，也可为 count×nodes 厚度温度剖面，供耦合
% 求解器在每个时间步从上一步完整剖面接续；不能只传表面温度。

%% 1. 网格与材料参数
if nargin<6, law=[]; end
backBoundary=wallBackBoundary(material);
if isfield(material,'layers') && ~isempty(material.layers)
    result=solveLayeredWallConduction(h,Taw,time,material,initial,law);
    return;
end
if backBoundary.active
    % 单层非绝热背面也使用同一守恒核，以完整 Newton 方程处理双面辐射。
    singleLayer=material;
    if ~isfield(singleLayer,'thermalNodes'), singleLayer.thermalNodes=12; end
    layered=struct('layers',{{singleLayer}},'contactResistance',[], ...
        'emissivity',material.emissivity,'radiationTemperature',material.radiationTemperature, ...
        'backBoundary',material.backBoundary);
    result=solveLayeredWallConduction(h,Taw,time,layered,initial,law);
    return;
end
if isfield(material,'propertyTable') && ~isempty(material.propertyTable)
    result=solveVariableMaterialConduction(h,Taw,time,material,initial,law);
    return;
end
nodes=12;
if isfield(material,'thermalNodes'), nodes=material.thermalNodes; end
validateattributes(nodes,{'numeric'},{'real','finite','scalar','integer','>=',3,'<=',100});
h=h(:); Taw=Taw(:); count=numel(h); time=time(:)'; nt=numel(time);
if numel(Taw)~=count
    error('RocketAeroHeat:Conduction1D:Size','h 与 Taw 尺寸不一致。');
end
dx=material.thickness/(nodes-1);
conductance=material.thermalConductivity/dx;
capacity=material.density*material.specificHeat*dx*ones(1,nodes);
capacity([1 end])=capacity([1 end])/2;
sigma=5.670374419e-8;
rad=material.emissivity*sigma;
Tr=material.radiationTemperature;
if isscalar(initial)
    validateattributes(initial,{'numeric'},{'real','finite','positive'});
    old=initial*ones(count,nodes);
elseif isequal(size(initial),[count nodes])
    validateattributes(initial,{'numeric'},{'real','finite','positive'});
    old=double(initial);
else
    error('RocketAeroHeat:Conduction1D:Initial','初温须为标量或面元数×厚度节点数的剖面。');
end
surface=nan(count,nt); surface(:,1)=old(:,1);
back=nan(count,nt); back(:,1)=old(:,end);
profile=old;

%% 2. 隐式有限体积步进与表面辐射 Newton 线性化
for step=2:nt
    dt=time(step)-time(step-1);
    guess=old(:,1);
    for iteration=1:25
        [currentH,derivativeH]=evaluateWallConvection(h,Taw,guess,law);
        flux=currentH.*(Taw-guess)-rad*(guess.^4-Tr^4);
        fluxDerivative=derivativeH.*(Taw-guess)-currentH-4*rad*guess.^3;
        diagonal=repmat(capacity/dt,count,1);
        diagonal(:,1)=diagonal(:,1)+conductance-fluxDerivative;
        diagonal(:,2:end-1)=diagonal(:,2:end-1)+2*conductance;
        diagonal(:,end)=diagonal(:,end)+conductance;
        rhs=old.*repmat(capacity/dt,count,1);
        rhs(:,1)=rhs(:,1)+flux-fluxDerivative.*guess;
        updated=solveTridiagonal(diagonal,rhs,conductance);
        if any(~isfinite(updated),'all') || any(updated(:)<=0)
            error('RocketAeroHeat:Conduction1D:Numerical', ...
                '一维导热求解得到非物理温度。');
        end
        difference=updated(:,1)-guess;
        guess=updated(:,1);
        if all(abs(difference)<=1e-8*max(1,guess)), break; end
    end
    if any(abs(difference)>1e-6*max(1,guess))
        error('RocketAeroHeat:Conduction1D:Convergence','表面辐射迭代未收敛。');
    end
    old=updated;
    surface(:,step)=updated(:,1);
    back(:,step)=updated(:,end);
end
profile=old;

%% 3. 保存方法和时间历程
result=struct;
result.method='一维有限体积导热 + 隐式 Euler + 温度相关对流/辐射 Newton 迭代';
result.scope='均匀常物性材料、固定外流、绝热背面；不含多层/烧蚀';
result.surfaceTemperature=surface;
result.backTemperature=back;
result.finalDepthTemperature=profile;
result.depth=linspace(0,material.thickness,nodes);
result.thermalNodes=nodes;
end

function answer=solveTridiagonal(diagonal,rhs,g)
% Thomas 法批量求解面元独立的常下/上对角三对角系统。
[count,n]=size(diagonal);
upper=zeros(count,n); forward=zeros(count,n);
upper(:,1)=-g./diagonal(:,1);
forward(:,1)=rhs(:,1)./diagonal(:,1);
for j=2:n
    denominator=diagonal(:,j)+g*upper(:,j-1);
    forward(:,j)=(rhs(:,j)+g*forward(:,j-1))./denominator;
    if j<n, upper(:,j)=-g./denominator; end
end
answer=zeros(count,n); answer(:,end)=forward(:,end);
for j=n-1:-1:1
    answer(:,j)=forward(:,j)-upper(:,j).*answer(:,j+1);
end
end
