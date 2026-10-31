function state = solveThermallyPerfectNormalShock(freestream)
%SOLVETHERMALLYPERFECTNORMALSHOCK 温度相关干空气的局部法向激波。
% 模块作用：以固定 N2/O2/Ar 组分 NASA7 焓与标准熵，在质量、动量、
% 总焓守恒和理想气体状态方程下反解法向激波后 T、p、rho、U、Mach。
% 再从同一焓/熵关系求等熵减速后的总温和总压损失。
% 可用于已验证球鼻驻点中心线，也可作二维切楔斜激波法向分量的
% 守恒内核；不求脱体距离、三维激波层、组分解离/电离或激波辐射。
% 任一状态超出 80~2500 K 时无解；80~150 K 用振动冻结近似。

%% 1. 检查来流状态，建立不能越界的输出契约
required={'Mach','gamma','R','pressure','temperature','density','speed'};
for index=1:numel(required)
    key=required{index};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:ThermalNormalShock:Input', ...
            '缺少 freestream.%s。',key);
    end
    validateattributes(freestream.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['freestream.' key]);
end
if freestream.Mach<=1
    error('RocketAeroHeat:ThermalNormalShock:Mach','法向激波要求 M>1。');
end
selectAirThermalModel(freestream,struct());
R=freestream.R; T1=freestream.temperature;
p1=freestream.pressure; rho1=freestream.density; U1=freestream.speed;
if abs(rho1*R*T1/p1-1)>1e-5 || ...
        abs(U1/(freestream.Mach*sqrt(freestream.gamma*R*T1))-1)>1e-5
    error('RocketAeroHeat:ThermalNormalShock:State', ...
        '来流密度或速度与输入状态方程不一致。');
end
state=struct('available',false, ...
    'method','固定组分干空气 NASA7/低温振动冻结焓熵 + 法向 Rankine-Hugoniot 跃迁', ...
    'scope',['单点法向跃迁；球鼻中心线或二维切楔法向分量；80~2500 K；' ...
        '固定组分、理想气体；不含解离、电离、激波辐射或三维激波层'], ...
    'thermalPropertyMode','thermallyPerfect', ...
    'status','尚未求解');
[cp1,h1,s1]=airThermalProperties(T1);
if ~isfinite(cp1) || cp1<=R
    state.status='来流温度超出干空气热物性适用域'; return;
end
gamma1=cp1/(cp1-R);
M1Thermal=U1/sqrt(gamma1*R*T1);
if M1Thermal<=1
    state.status='温度相关比热下的来流不满足超声速条件'; return;
end
T0=airRecoveryTemperature(T1,U1,1);
if ~isfinite(T0)
    state.status='总焓对应温度超过干空气热物性上限'; return;
end
[~,h0,s0]=airThermalProperties(T0);

%% 2. 以密度比求非平凡压缩根，弱激波也不能丢掉
massFlux=rho1*U1;
% 旧的温度-慢速度根只覆盖强压缩分支：M_n 约 1~1.2 时，真实的
% 亚声速激波后速度仍在二次方程的快根上。改以 r=rho2/rho1>1：
% U2=U1/r；动量方程给 p2，再由状态方程求 T2。r=1 是平凡
% 无激波根，超声速来流在其右侧能量残差先正后负，二分取另一根。
if M1Thermal-1<1e-8
    state.status='法向 Mach 数过于接近 1，激波强度低于求解精度'; return;
end
lower=1+min(1e-4,0.25*(M1Thermal-1));
upper=1e4;
if compressionEnergy(lower,T1,U1,R,h0)<=0 || ...
        compressionEnergy(upper,T1,U1,R,h0)>=0
    state.status='法向激波压缩根没有有效括区'; return;
end
for iteration=1:72
    trial=0.5*(lower+upper);
    if compressionEnergy(trial,T1,U1,R,h0)>0
        lower=trial;
    else
        upper=trial;
    end
end
densityRatio=0.5*(lower+upper);
U2=U1/densityRatio;
p2=p1+rho1*U1*(U1-U2);
rho2=rho1*densityRatio;
T2=p2/(rho2*R);
if T2<80 || T2>2500
    state.status='激波后温度超出干空气热物性适用域'; return;
end
[cp2,h2,s2]=airThermalProperties(T2);
gamma2=cp2/(cp2-R);
M2=U2/sqrt(gamma2*R*T2);
if ~isfinite(M2) || M2>=1 || p2<=p1 || T2<=T1
    state.status='未得到物理压缩且亚声速的激波后状态'; return;
end

%% 3. 由同一标准熵积分计算前后总压和守恒残差
p01=p1*exp((s0-s1)/R);
p02=p2*exp((s0-s2)/R);
entropyIncrease=(s2-s1)-R*log(p2/p1);
massResidual=(rho2*U2-massFlux)/massFlux;
momentum1=p1+rho1*U1^2;
momentumResidual=(p2+rho2*U2^2-momentum1)/momentum1;
energyResidual=(h2+0.5*U2^2-h0)/max(1,abs(h0));
if ~isfinite(p02) || p02<=0 || p02>=p01 || ...
        entropyIncrease<=0 || ...
        max(abs([massResidual,momentumResidual,energyResidual]))>1e-8
    state.status='质量、动量、总焓或熵增校核未通过'; return;
end
state.available=true;
state.status='局部温度相关法向激波后状态已计算；激波几何仍未求解';
state.upstreamMach=freestream.Mach;
state.upstreamMachThermal=M1Thermal;
state.postShockMach=M2;
state.postShockPressure=p2;
state.postShockTemperature=T2;
state.postShockDensity=rho2;
state.postShockSpeed=U2;
state.postShockTotalPressure=p02;
state.postShockTotalTemperature=T0;
state.stagnationPressure=p02;
state.stagnationTemperature=T0;
state.upstreamTotalPressureThermal=p01;
state.pressureRatio=p2/p1;
state.densityRatio=rho2/rho1;
state.temperatureRatio=T2/T1;
state.totalPressureRatio=p02/p01;
state.entropyIncrease=entropyIncrease;
state.shockStandOffDistance=NaN;
state.shockSurfacePosition=[NaN NaN NaN];
state.massResidual=massResidual;
state.momentumResidual=momentumResidual;
state.energyResidual=energyResidual;
end

function residual=compressionEnergy(r,T1,U1,R,h0)
% 质量/动量/状态方程消元后的总焓差；越界温度只用符号括区，不外推。
T2=T1/r+(U1^2/R)*(r-1)/r^2;
if T2>2500
    residual=1; return;
elseif T2<80
    residual=-1; return;
end
[~,h2]=airThermalProperties(T2);
residual=h2+0.5*(U1/r)^2-h0;
end
