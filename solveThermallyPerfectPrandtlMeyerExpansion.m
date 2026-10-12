function state=solveThermallyPerfectPrandtlMeyerExpansion(upstream,turnDeg)
%SOLVETHERMALLYPERFECTPRANDTLMEYEREXPANSION 干空气二维等熵膨胀扇。
% 模块作用：固定组分 NASA7 焓/熵给出局部音速与等熵压力；积分
% d(theta)=sqrt(M^2-1)*dU/U，再按指定凸转角反求下游温度和状态。
% 仅适用于局部二维超声速简单波，不能表示三维波系相互作用或分离。
% 静温/总温必须处于 80~2500 K；80~150 K 用振动冻结近似。

%% 1. 校验上游流动，建立与旧版一致的结构体接口
required={'Mach','gamma','R','pressure','temperature','density','speed'};
if ~isstruct(upstream) || ~isscalar(upstream)
    error('RocketAeroHeat:ThermalExpansion:Input','上游状态须为标量结构体。');
end
for k=1:numel(required)
    key=required{k};
    if ~isfield(upstream,key)
        error('RocketAeroHeat:ThermalExpansion:Input','缺少 upstream.%s。',key);
    end
    validateattributes(upstream.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['upstream.' key]);
end
validateattributes(turnDeg,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<',90},mfilename,'turnDeg');
selectAirThermalModel(upstream,struct());
R=upstream.R; T1=upstream.temperature;
p1=upstream.pressure; U1=upstream.speed;
if abs(upstream.density*R*T1/p1-1)>1e-5 || ...
        abs(U1/(upstream.Mach*sqrt(upstream.gamma*R*T1))-1)>1e-5
    error('RocketAeroHeat:ThermalExpansion:State', ...
        '上游密度或速度与输入状态方程不一致。');
end
state=struct('available',false, ...
    'method','固定组分干空气 NASA7/低温振动冻结焓熵二维 Prandtl-Meyer 膨胀', ...
    'scope',['局部二维超声速凸转角简单波；80~2500 K；' ...
        '不含三维波系相互作用、分离、尾迹或化学反应'], ...
    'thermalPropertyMode','thermallyPerfect','turnDeg',turnDeg, ...
    'status','尚未求解');
fields={'maximumAdditionalTurnDeg','upstreamMachThermal','Mach', ...
    'pressure','temperature','density','speed','speedOfSound', ...
    'totalPressure','totalTemperature','pressureRatio', ...
    'temperatureRatio','entropyResidual','energyResidual'};
for k=1:numel(fields), state.(fields{k})=NaN; end
[cp1,h1,s1]=airThermalProperties(T1);
T0=airRecoveryTemperature(T1,U1,1);
if ~isfinite(cp1) || cp1<=R || ~isfinite(T0)
    state.status='上游静温或总焓超出干空气热物性适用域'; return;
end
[~,h0,s0]=airThermalProperties(T0);
a1=sqrt(cp1/(cp1-R)*R*T1);
M1=U1/a1;
if M1<=1 || upstream.Mach<=1
    state.status='温度相关比热下的上游流动不是超声速'; return;
end
p0=p1*exp((s0-s1)/R);
state.upstreamMachThermal=M1;
state.totalPressure=p0;
state.totalTemperature=T0;

%% 2. 在当前上游总焓下建立膨胀角—温度曲线
% 同一上游温度和速度的多个面元共用积分曲线；压力只用于最后
% 的等熵比例，故不属于缓存键。4,001 节点使角度插值误差很小。
persistent cachedKey cachedCurve
key=[T1,U1,R];
if isempty(cachedKey) || ~isequal(cachedKey,key)
    grid=linspace(T1,80,4001)';
    [cp,h,s]=airThermalProperties(grid);
    gamma=cp./(cp-R);
    velocity=sqrt(max(0,2*(h0-h)));
    sound=sqrt(gamma.*R.*grid);
    mach=velocity./sound;
    integrand=sqrt(max(0,mach.^2-1)).*cp./velocity.^2;
    angle=cumtrapz(T1-grid,integrand);
    cachedKey=key;
    cachedCurve=struct('temperature',grid,'angle',angle, ...
        'maximum',angle(end));
end
curve=cachedCurve;
state.maximumAdditionalTurnDeg=rad2deg(curve.maximum);
target=deg2rad(turnDeg);

%% 3. 零转角保持上游状态；超出温度下限时不给虚假膨胀解
if turnDeg==0
    state.available=true;
    state.status='零转角：保持上游等熵状态';
    state.Mach=M1; state.pressure=p1; state.temperature=T1;
    state.density=upstream.density; state.speed=U1;
    state.speedOfSound=a1;
    state.pressureRatio=1; state.temperatureRatio=1;
    state.entropyResidual=0; state.energyResidual=0;
    return;
end
if target>curve.maximum+1e-10
    state.status='指定膨胀角要求静温低于 80 K，热物性域内无解'; return;
end
T2=interp1(curve.angle,curve.temperature,target,'pchip');
[cp2,h2,s2]=airThermalProperties(T2);
p2=p1*exp((s2-s1)/R);
rho2=p2/(R*T2);
U2=sqrt(max(0,2*(h0-h2)));
a2=sqrt(cp2/(cp2-R)*R*T2);
M2=U2/a2;
entropyResidual=(s2-s1)-R*log(p2/p1);
energyResidual=(h2+0.5*U2^2-h0)/max(1,abs(h0));
if ~isfinite(T2) || T2<80 || T2>T1 || ...
        p2<=0 || p2>=p1 || M2<=M1 || ...
        abs(entropyResidual)>1e-7 || abs(energyResidual)>1e-9
    state.status='膨胀后的等熵、能量或超声速校核未通过'; return;
end
state.available=true;
state.status='局部二维温度相关等熵膨胀解';
state.Mach=M2; state.pressure=p2; state.temperature=T2;
state.density=rho2; state.speed=U2; state.speedOfSound=a2;
state.pressureRatio=p2/p1;
state.temperatureRatio=T2/T1;
state.entropyResidual=entropyResidual;
state.energyResidual=energyResidual;
end
