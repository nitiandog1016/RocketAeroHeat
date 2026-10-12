function state = solvePrandtlMeyerExpansion(freestream,turnDeg)
%SOLVEPRANDTLMEYEREXPANSION 二维定比热完全气体等熵膨胀扇。
% 输入上游静态/总状态和正的背风转折角 [deg]，返回膨胀后的均匀状态。
% 仅适用于超声速、局部凸转角；不描述真实三维尾迹或分离流。

%% 1. 验证输入及可达到的最大膨胀角
validateattributes(turnDeg,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<',90},mfilename,'turnDeg');
for name={'Mach','gamma','R','pressure','temperature','density', ...
        'totalPressure','totalTemperature'}
    key=name{1};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:Expansion:Input','缺少 freestream.%s。',key);
    end
    validateattributes(freestream.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['freestream.' key]);
end
M1=freestream.Mach; g=freestream.gamma;
if M1<=1 || g<=1
    error('RocketAeroHeat:Expansion:Domain','膨胀扇要求 M>1 且 gamma>1。');
end
nu1=prandtlMeyer(M1,g);
nuMaximum=pi/2*(sqrt((g+1)/(g-1))-1);
target=nu1+deg2rad(turnDeg);
state=struct;
state.method='二维 Prandtl-Meyer 等熵膨胀扇';
state.scope='局部凸转角、超声速完全气体；未描述三维肩部干扰和分离';
state.turnDeg=turnDeg; state.maximumAdditionalTurnDeg=rad2deg(nuMaximum-nu1);
state.available=target<nuMaximum;
for name={'Mach','pressure','temperature','density','speed', ...
        'totalPressure','totalTemperature','pressureRatio','temperatureRatio'}
    state.(name{1})=NaN;
end
if ~state.available
    state.status='膨胀角超出有限马赫数的 Prandtl-Meyer 范围';
    return;
end

%% 2. 单调反求下游马赫数；零转折精确返回上游状态
if turnDeg==0
    M2=M1;
else
    upper=max(2*M1,2);
    while prandtlMeyer(upper,g)<target
        upper=2*upper;
        if upper>1e8
            state.available=false;
            state.status='膨胀角接近渐近极限，有限精度下无法求解';
            return;
        end
    end
    M2=fzero(@(m)prandtlMeyer(m,g)-target,[M1,upper], ...
        optimset('TolX',1e-12,'Display','off'));
end

%% 3. 等熵关系恢复下游状态，总压和总温保持不变
factor=(1+0.5*(g-1)*M1^2)/(1+0.5*(g-1)*M2^2);
T2=freestream.temperature*factor;
p2=freestream.pressure*factor^(g/(g-1));
rho2=p2/(freestream.R*T2);
state.Mach=M2; state.pressure=p2; state.temperature=T2;
state.density=rho2; state.speed=M2*sqrt(g*freestream.R*T2);
state.totalPressure=freestream.totalPressure;
state.totalTemperature=freestream.totalTemperature;
state.pressureRatio=p2/freestream.pressure;
state.temperatureRatio=T2/freestream.temperature;
state.status='局部等熵膨胀解';
end

function nu=prandtlMeyer(M,gamma)
root=sqrt(M^2-1);
nu=sqrt((gamma+1)/(gamma-1))* ...
    atan(sqrt((gamma-1)/(gamma+1))*root)-atan(root);
end
