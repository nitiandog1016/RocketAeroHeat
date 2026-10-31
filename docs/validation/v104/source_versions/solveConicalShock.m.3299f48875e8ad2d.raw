function state = solveConicalShock(freestream,coneHalfAngleDeg)
%SOLVECONICALSHOCK 零攻角、直圆锥的 Taylor-Maccoll 附体锥形激波弱解。
% 输入为自由来流和几何半锥角 [deg]。输出为锥面边界层外缘状态。
% 与二维直楔斜激波不同，本模块积分轴对称无粘锥流常微分方程。
% 不适用于曲面鼻锥、非零攻角、肩部干扰、脱体激波或真实气体。

%% 1. 输入和方法说明
validateattributes(coneHalfAngleDeg,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<',90},mfilename,'coneHalfAngleDeg');
for name={'Mach','gamma','R','pressure','temperature','density','speed', ...
        'totalPressure','totalTemperature'}
    key=name{1};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:ConicalShock:Input','缺少 freestream.%s。',key);
    end
    validateattributes(freestream.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['freestream.' key]);
end
M1=freestream.Mach; gamma=freestream.gamma;
if M1<=1 || gamma<=1
    error('RocketAeroHeat:ConicalShock:Domain','Taylor-Maccoll 附体解要求 M>1、gamma>1。');
end
state=struct;
state.method='Taylor-Maccoll 轴对称无粘锥流 + 正激波法向跃迁';
state.scope='直圆锥、零攻角、定比热完全气体、附体弱解';
state.coneHalfAngleDeg=coneHalfAngleDeg;
state.upstreamMach=M1;
state.attached=false;
for name={'shockAngleDeg','Mach','pressure','temperature','density','speed', ...
        'totalPressure','totalTemperature','pressureRatio','totalPressureRatio', ...
        'radialVelocityAtWall','polarVelocityAtWall'}
    state.(name{1})=NaN;
end
if coneHalfAngleDeg==0
    state.attached=true; state.status='零半锥角：自由来流，无激波';
    state.shockAngleDeg=asind(1/M1);
    state.Mach=M1; state.pressure=freestream.pressure;
    state.temperature=freestream.temperature; state.density=freestream.density;
    state.speed=freestream.speed; state.totalPressure=freestream.totalPressure;
    state.totalTemperature=freestream.totalTemperature;
    state.pressureRatio=1; state.totalPressureRatio=1;
    state.radialVelocityAtWall=freestream.speed;
    state.polarVelocityAtWall=0;
    return;
end

%% 2. 扫描激波角并寻找最小的附体根（弱解）
thetaC=deg2rad(coneHalfAngleDeg); machAngle=asin(1/M1);
lower=max(machAngle,thetaC)+1e-7;
upper=pi/2-1e-5;
if lower>=upper
    state.status='半锥角超出数值附体范围'; return;
end
angles=linspace(lower,upper,80);
residual=nan(size(angles));
for i=1:numel(angles)
    residual(i)=wallPolarVelocity(angles(i),thetaC,M1,gamma);
end
pair=find(isfinite(residual(1:end-1)) & isfinite(residual(2:end)) & ...
    residual(1:end-1).*residual(2:end)<=0,1,'first');
if isempty(pair)
    state.status='未找到附体锥形激波弱解；可能已脱体'; return;
end
beta=fzero(@(b)wallPolarVelocity(b,thetaC,M1,gamma), ...
    [angles(pair),angles(pair+1)],optimset('TolX',1e-11,'Display','off'));
[wall,normalRatio]=integrateCone(beta,thetaC,M1,gamma);
if any(~isfinite(wall)) || abs(wall(2))>2e-5
    state.status='锥面边界条件未收敛'; return;
end

%% 3. 由激波后等熵锥流状态恢复锥面热力参数
Mn1=M1*sin(beta); nsq=Mn1^2;
p2Ratio=1+2*gamma/(gamma+1)*(nsq-1);
rho2Ratio=normalRatio;
Mn2=sqrt((1+0.5*(gamma-1)*nsq)/(gamma*nsq-0.5*(gamma-1)));
deflection=atan(2*(nsq-1)/tan(beta)/ ...
    (M1^2*(gamma+cos(2*beta))+2));
M2=Mn2/sin(beta-deflection);
p02=freestream.pressure*p2Ratio* ...
    (1+0.5*(gamma-1)*M2^2)^(gamma/(gamma-1));
u2=wall(1)^2+wall(2)^2;
if u2>=1 || u2<=0
    state.status='积分结果超出气体总焓范围'; return;
end
Tw=freestream.totalTemperature*(1-u2);
Mw=sqrt(2*u2/((gamma-1)*(1-u2)));
pw=p02/(1+0.5*(gamma-1)*Mw^2)^(gamma/(gamma-1));
rw=pw/(freestream.R*Tw);
speed=Mw*sqrt(gamma*freestream.R*Tw);
state.attached=true; state.status='附体弱锥形激波；锥面状态已收敛';
state.shockAngleDeg=rad2deg(beta);
state.Mach=Mw; state.pressure=pw; state.temperature=Tw;
state.density=rw; state.speed=speed;
state.totalPressure=p02; state.totalTemperature=freestream.totalTemperature;
state.pressureRatio=pw/freestream.pressure;
state.totalPressureRatio=p02/freestream.totalPressure;
normalizer=sqrt(2*gamma*freestream.R*freestream.totalTemperature/(gamma-1));
state.radialVelocityAtWall=wall(1)*normalizer;
state.polarVelocityAtWall=wall(2)*normalizer;
end

function polar=wallPolarVelocity(beta,thetaC,M,gamma)
[wall,~]=integrateCone(beta,thetaC,M,gamma);
polar=wall(2);
end

function [velocity,densityRatio]=integrateCone(beta,thetaC,M,gamma)
% 四阶 Runge-Kutta 从锥形激波面向锥面反向积分，遇奇异点返回 NaN。
Mn=M*sin(beta);
densityRatio=(gamma+1)*Mn^2/((gamma-1)*Mn^2+2);
normalizedSpeed=sqrt((gamma-1)*M^2/(2+(gamma-1)*M^2));
velocity=[normalizedSpeed*cos(beta);-normalizedSpeed*sin(beta)/densityRatio];
steps=max(100,ceil(rad2deg(beta-thetaC)*8));
h=(thetaC-beta)/steps;
for j=1:steps
    angle=beta+(j-1)*h;
    k1=coneDerivative(angle,velocity,gamma);
    k2=coneDerivative(angle+h/2,velocity+h*k1/2,gamma);
    k3=coneDerivative(angle+h/2,velocity+h*k2/2,gamma);
    k4=coneDerivative(angle+h,velocity+h*k3,gamma);
    if any(~isfinite([k1;k2;k3;k4]))
        velocity=[NaN;NaN]; return;
    end
    velocity=velocity+h*(k1+2*k2+2*k3+k4)/6;
    if ~all(isfinite(velocity)) || sum(velocity.^2)>=1
        velocity=[NaN;NaN]; return;
    end
end
end

function derivative=coneDerivative(theta,v,gamma)
radial=v(1); polar=v(2);
thermal=(gamma-1)/2*(1-radial^2-polar^2);
denominator=thermal-polar^2;
if abs(denominator)<1e-10 || thermal<=0 || theta<=0
    derivative=[NaN;NaN]; return;
end
derivative=[polar; (polar^2*radial- ...
    thermal*(2*radial+polar/tan(theta)))/denominator];
end
