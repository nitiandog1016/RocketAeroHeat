function shock=solveThermallyPerfectObliqueShock(upstream,deflectionDeg)
%SOLVETHERMALLYPERFECTOBLIQUESHOCK 固定组分干空气的二维附体弱斜激波。
% 模块作用：由激波法向 Rankine-Hugoniot 守恒、切向速度连续及 NASA7
% 干空气焓/熵求激波角和下游状态。总压用同一等熵熵差反算。
% 仅适用于局部二维附体直楔；不求三维激波面、锥流、脱体波或化学反应。
% 全部静温与总温须在 80~2500 K，超出时返回无解。

%% 1. 校验输入并建立明确的无解输出
required={'Mach','gamma','R','pressure','temperature','density','speed'};
if ~isstruct(upstream) || ~isscalar(upstream)
    error('RocketAeroHeat:ThermalObliqueShock:Input','上游状态须为标量结构体。');
end
for k=1:numel(required)
    key=required{k};
    if ~isfield(upstream,key)
        error('RocketAeroHeat:ThermalObliqueShock:Input','缺少 upstream.%s。',key);
    end
    validateattributes(upstream.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['upstream.' key]);
end
validateattributes(deflectionDeg,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<',90},mfilename,'deflectionDeg');
selectAirThermalModel(upstream,struct());
R=upstream.R; T1=upstream.temperature;
p1=upstream.pressure; rho1=upstream.density; U1=upstream.speed;
if abs(rho1*R*T1/p1-1)>1e-5 || ...
        abs(U1/(upstream.Mach*sqrt(upstream.gamma*R*T1))-1)>1e-5
    error('RocketAeroHeat:ThermalObliqueShock:State', ...
        '上游密度或速度与输入状态方程不一致。');
end
shock=struct('available',false,'attached',false,'branch','weak', ...
    'method','固定组分干空气 NASA7/低温振动冻结焓熵 + 二维附体弱斜激波', ...
    'scope',['局部二维直楔；80~2500 K；不含三维激波层、锥流、' ...
        '脱体激波、解离、电离或激波辐射'], ...
    'thermalPropertyMode','thermallyPerfect', ...
    'deflectionDeg',deflectionDeg,'status','尚未求解');
blank={'shockAngleDeg','machAngleDeg','peakShockAngleDeg', ...
    'maximumAttachedDeflectionDeg','upstreamMachThermal', ...
    'normalMachBefore','normalMachAfter','Mach','pressure', ...
    'temperature','density','speed','speedOfSound','totalPressure', ...
    'totalTemperature','upstreamTotalPressureThermal','pressureRatio', ...
    'temperatureRatio','densityRatio','totalPressureRatio', ...
    'entropyIncrease','massResidual','momentumResidual','energyResidual', ...
    'tangentialSpeedBefore','tangentialSpeedAfter'};
for k=1:numel(blank), shock.(blank{k})=NaN; end
[cp1,h1,s1]=airThermalProperties(T1);
T0=airRecoveryTemperature(T1,U1,1);
if ~isfinite(cp1) || cp1<=R || ~isfinite(T0)
    shock.status='上游静温或总焓超出干空气热物性适用域'; return;
end
[~,~,s0]=airThermalProperties(T0);
gamma1=cp1/(cp1-R);
a1=sqrt(gamma1*R*T1);
M1=U1/a1;
if M1<=1 || upstream.Mach<=1
    shock.status='温度相关比热下的上游流动不是超声速'; return;
end
p01=p1*exp((s0-s1)/R);
% 正激波核同时核对输入的定 gamma Mach；下界取两个声速中较大者。
legacyA=sqrt(upstream.gamma*R*T1);
betaMin=asin(max(a1,legacyA)/U1);
if ~isreal(betaMin) || betaMin>=pi/2
    shock.status='上游法向速度没有超声速附体范围'; return;
end
shock.upstreamMachThermal=M1;
shock.upstreamTotalPressureThermal=p01;
shock.totalTemperature=T0;
shock.machAngleDeg=rad2deg(betaMin);

%% 2. 求激波极线的最大附体转角；同一上游状态复用极限角
persistent cachedKey cachedPeak cachedMaximum
key=[T1,U1,R,upstream.gamma];
if isempty(cachedKey) || ~isequal(cachedKey,key)
    [peak,negativeMaximum]=fminbnd( ...
        @(beta)-turnAtBeta(beta,upstream,betaMin), ...
        betaMin,pi/2,optimset('TolX',1e-9,'Display','off'));
    cachedKey=key;
    cachedPeak=peak;
    cachedMaximum=-negativeMaximum;
end
peak=cachedPeak; maximum=cachedMaximum;
shock.peakShockAngleDeg=rad2deg(peak);
shock.maximumAttachedDeflectionDeg=rad2deg(maximum);
theta=deg2rad(deflectionDeg);
if theta>maximum+1e-9
    shock.status='脱体：二维温度相关附体弱斜激波解不存在'; return;
end
shock.attached=true;

%% 3. 零转角保持上游状态；非零转角反求弱支激波角
if theta==0
    shock.available=true;
    shock.status='零转角：自由来流，无激波';
    shock.shockAngleDeg=rad2deg(betaMin);
    shock.normalMachBefore=1; shock.normalMachAfter=1;
    shock.Mach=M1; shock.pressure=p1; shock.temperature=T1;
    shock.density=rho1; shock.speed=U1; shock.speedOfSound=a1;
    shock.totalPressure=p01;
    shock.pressureRatio=1; shock.temperatureRatio=1;
    shock.densityRatio=1; shock.totalPressureRatio=1;
    shock.entropyIncrease=0;
    shock.massResidual=0; shock.momentumResidual=0;
    shock.energyResidual=0;
    shock.tangentialSpeedBefore=U1*cos(betaMin);
    shock.tangentialSpeedAfter=shock.tangentialSpeedBefore;
    return;
end
if theta>=maximum-1e-9
    beta=peak;
else
    beta=fzero(@(angle)turnAtBeta(angle,upstream,betaMin)-theta, ...
        [betaMin,peak],optimset('TolX',1e-10,'Display','off'));
end
normal=normalJump(beta,upstream,betaMin);
if ~normal.available
    shock.attached=false;
    shock.status='激波后温度超出热物性域或守恒方程无解'; return;
end

%% 4. 保持切向速度连续，恢复全速度与热力总状态
Un1=U1*sin(beta); Ut1=U1*cos(beta);
Un2=normal.postShockSpeed;
U2=hypot(Un2,Ut1);
T2=normal.postShockTemperature; p2=normal.postShockPressure;
rho2=normal.postShockDensity;
[cp2,h2,s2]=airThermalProperties(T2);
gamma2=cp2/(cp2-R); a2=sqrt(gamma2*R*T2);
p02=p2*exp((s0-s2)/R);
entropyIncrease=(s2-s1)-R*log(p2/p1);
thetaCalculated=beta-atan2(Un2,Ut1);
massResidual=(rho2*Un2-rho1*Un1)/(rho1*Un1);
momentumResidual=(p2+rho2*Un2^2-p1-rho1*Un1^2)/ ...
    (p1+rho1*Un1^2);
energyResidual=(h2+0.5*U2^2-h1-0.5*U1^2)/ ...
    max(1,abs(h1+0.5*U1^2));
if ~isfinite(p02) || p02<=0 || p02>=p01 || ...
        entropyIncrease<=0 || abs(thetaCalculated-theta)>1e-7 || ...
        max(abs([massResidual,momentumResidual,energyResidual]))>1e-8
    shock.attached=false;
    shock.status='斜激波守恒、转角或熵增校核未通过'; return;
end
shock.available=true;
shock.status='二维局部附体弱斜激波后热力状态已计算';
shock.shockAngleDeg=rad2deg(beta);
shock.normalMachBefore=Un1/a1;
shock.normalMachAfter=Un2/a2;
shock.Mach=U2/a2;
shock.pressure=p2; shock.temperature=T2; shock.density=rho2;
shock.speed=U2; shock.speedOfSound=a2;
shock.totalPressure=p02;
shock.pressureRatio=p2/p1;
shock.temperatureRatio=T2/T1;
shock.densityRatio=rho2/rho1;
shock.totalPressureRatio=p02/p01;
shock.entropyIncrease=entropyIncrease;
shock.massResidual=massResidual;
shock.momentumResidual=momentumResidual;
shock.energyResidual=energyResidual;
shock.tangentialSpeedBefore=Ut1;
shock.tangentialSpeedAfter=Ut1;
end

function theta=turnAtBeta(beta,upstream,betaMin)
% 激波极线上的转角；Mach 角或无热力解处返回零。
if beta<=betaMin+1e-12
    theta=0; return;
end
state=normalJump(beta,upstream,betaMin);
if ~state.available
    theta=0; return;
end
theta=beta-atan2(state.postShockSpeed,upstream.speed*cos(beta));
end

function state=normalJump(beta,upstream,betaMin)
% 只向现有法向守恒核传入速度法向分量；全流总压在主函数另算。
state=struct('available',false);
if beta<=betaMin+1e-12, return; end
stateInput=upstream;
stateInput.speed=upstream.speed*sin(beta);
stateInput.Mach=stateInput.speed/ ...
    sqrt(upstream.gamma*upstream.R*upstream.temperature);
if stateInput.Mach<=1, return; end
state=solveThermallyPerfectNormalShock(stateInput);
end
