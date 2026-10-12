function shock = solveObliqueShock(freestream,deflectionDeg)
%SOLVEOBLIQUESHOCK 定比热完全气体的二维附体斜激波（弱解）。
% 输入：freestream 为 calcFreestream 的输出；deflectionDeg 为流动转折角 [deg]。
% 输出：shock 的各物理量均为激波后、激波外侧的均匀状态。
% 本模块求解二维直楔的 theta-beta-M 关系；不将其冒充三维锥流或整箭流场。
% 当转折角超过附体极限时，status='detached'，激波角和后状态为 NaN。
% theta=0 精确返回自由来流，不人为制造无限弱激波的熵增。

%% 1. 验证来流和角度
keys={'Mach','gamma','R','pressure','temperature','density','speed', ...
    'totalPressure','totalTemperature'};
if ~isstruct(freestream) || ~isscalar(freestream)
    error('RocketAeroHeat:ObliqueShock:Input','freestream 必须为标量结构体。');
end
for i=1:numel(keys)
    name=keys{i};
    if ~isfield(freestream,name)
        error('RocketAeroHeat:ObliqueShock:Input','缺少 freestream.%s。',name);
    end
    validateattributes(freestream.(name),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['freestream.' name]);
end
validateattributes(deflectionDeg,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<',90},mfilename,'deflectionDeg');
M1=double(freestream.Mach); gamma=double(freestream.gamma);
if M1<=1 || gamma<=1
    error('RocketAeroHeat:ObliqueShock:Domain','附体斜激波要求 M>1、gamma>1。');
end

%% 2. 求附体极限，区分弱解与脱体情形
machAngle=asin(1/M1); normalAngle=pi/2;
% 在 Mach 角与正激波角之间寻找最大可附体转折角。
[peakBeta,negativePeak]=fminbnd(@(beta)-turnAngle(beta,M1,gamma), ...
    machAngle,normalAngle,optimset('TolX',1e-12,'Display','off'));
maximumAngle=-negativePeak;
theta=deg2rad(deflectionDeg);
shock=struct;
shock.method='二维定比热完全气体 theta-beta-M 关系（附体弱斜激波）';
shock.scope='局部二维直楔；不是轴对称锥流、三维激波形状或驻点解';
shock.upstreamMach=M1; shock.gamma=gamma;
shock.deflectionDeg=deflectionDeg;
shock.machAngleDeg=rad2deg(machAngle);
shock.maximumAttachedDeflectionDeg=rad2deg(maximumAngle);
shock.peakShockAngleDeg=rad2deg(peakBeta);
shock.attached=theta<=maximumAngle+1e-11;
shock.branch='weak';
% 先填 NaN，脱体时避免留有看似可信的下游状态。
names={'shockAngleDeg','normalMachBefore','normalMachAfter','Mach', ...
    'pressure','temperature','density','speed','speedOfSound', ...
    'totalPressure','totalTemperature','pressureRatio','temperatureRatio', ...
    'densityRatio','totalPressureRatio','entropyIncrease'};
for i=1:numel(names), shock.(names{i})=NaN; end
if ~shock.attached
    shock.status='脱体：二维附体斜激波弱解不存在；需要脱体激波/外流求解';
    return;
end

%% 3. 在弱分支反求激波角，转折零时保持自由来流
if theta==0
    beta=machAngle;
    shock.status='零转折：自由来流，无激波';
    shock.shockAngleDeg=rad2deg(beta);
    shock.normalMachBefore=1; shock.normalMachAfter=1;
    shock.Mach=M1; shock.pressure=freestream.pressure;
    shock.temperature=freestream.temperature;
    shock.density=freestream.density;
    shock.speed=freestream.speed;
    shock.speedOfSound=sqrt(gamma*freestream.R*shock.temperature);
    shock.totalPressure=freestream.totalPressure;
    shock.totalTemperature=freestream.totalTemperature;
    shock.pressureRatio=1; shock.temperatureRatio=1;
    shock.densityRatio=1; shock.totalPressureRatio=1;
    shock.entropyIncrease=0;
    return;
end
if theta>=maximumAngle-1e-11
    beta=peakBeta;
else
    beta=fzero(@(b)turnAngle(b,M1,gamma)-theta, ...
        [machAngle,peakBeta],optimset('TolX',1e-12,'Display','off'));
end

%% 4. 在激波法向使用 Rankine-Hugoniot 跃迁，再恢复下游流向
Mn1=M1*sin(beta); Mn1sq=Mn1^2;
pressureRatio=1+2*gamma/(gamma+1)*(Mn1sq-1);
densityRatio=(gamma+1)*Mn1sq/((gamma-1)*Mn1sq+2);
temperatureRatio=pressureRatio/densityRatio;
Mn2sq=(1+0.5*(gamma-1)*Mn1sq)/(gamma*Mn1sq-0.5*(gamma-1));
Mn2=sqrt(Mn2sq); M2=Mn2/sin(beta-theta);
p2=freestream.pressure*pressureRatio;
T2=freestream.temperature*temperatureRatio;
rho2=freestream.density*densityRatio;
a2=sqrt(gamma*freestream.R*T2);
p02=p2*(1+0.5*(gamma-1)*M2^2)^(gamma/(gamma-1));
T02=T2*(1+0.5*(gamma-1)*M2^2);

%% 5. 保存算法、适用范围和可直接校核的跃迁比
shock.shockAngleDeg=rad2deg(beta);
shock.normalMachBefore=Mn1; shock.normalMachAfter=Mn2;
shock.Mach=M2; shock.pressure=p2; shock.temperature=T2;
shock.density=rho2; shock.speed=M2*a2; shock.speedOfSound=a2;
shock.totalPressure=p02; shock.totalTemperature=T02;
shock.pressureRatio=pressureRatio; shock.temperatureRatio=temperatureRatio;
shock.densityRatio=densityRatio;
shock.totalPressureRatio=p02/freestream.totalPressure;
shock.entropyIncrease=freestream.R*log(freestream.totalPressure/p02);
shock.status='附体弱解；下游状态有效，但仅适用于局部二维直楔';
end

function theta=turnAngle(beta,M,gamma)
% theta-beta-M 关系：用 atan2 保持弱激波区的数值稳定。
normal=M^2*sin(beta)^2;
numerator=2*(normal-1)/tan(beta);
denominator=M^2*(gamma+cos(2*beta))+2;
theta=atan2(numerator,denominator);
end
