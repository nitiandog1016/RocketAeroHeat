function state = solveNormalShock(freestream)
%SOLVENORMALSHOCK 完全气体法向激波的局部 Rankine-Hugoniot 跃迁。
% 模块作用：给球形钝鼻轴向驻点前的激波中心线提供热力状态。
% 计算激波后 p、T、rho、U、M、总压损失和等熵减速后的驻点状态。
% 只确定局部法向跃迁；脱体激波距离、三维形状及鼻部侧向流场仍无解。

%% 1. 检查自由来流与完全气体适用域
required={'Mach','gamma','R','pressure','temperature','density', ...
    'speed','totalPressure','totalTemperature'};
for j=1:numel(required)
    key=required{j};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:NormalShock:Input', ...
            '缺少 freestream.%s。',key);
    end
    validateattributes(freestream.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['freestream.' key]);
end
M1=freestream.Mach; gamma=freestream.gamma;
if M1<=1 || gamma<=1
    error('RocketAeroHeat:NormalShock:Domain', ...
        '法向激波要求 M>1 且 gamma>1。');
end
expectedDensity=freestream.pressure/(freestream.R*freestream.temperature);
expectedSpeed=M1*sqrt(gamma*freestream.R*freestream.temperature);
if abs(freestream.density/expectedDensity-1)>1e-5 || ...
        abs(freestream.speed/expectedSpeed-1)>1e-5
    error('RocketAeroHeat:NormalShock:State', ...
        '来流密度、速度与完全气体状态方程不一致。');
end

%% 2. 法向动量、质量、能量跃迁
M1s=M1^2;
pressureRatio=1+2*gamma/(gamma+1)*(M1s-1);
densityRatio=(gamma+1)*M1s/((gamma-1)*M1s+2);
temperatureRatio=pressureRatio/densityRatio;
M2=sqrt((1+0.5*(gamma-1)*M1s)/ ...
    (gamma*M1s-0.5*(gamma-1)));
p2=freestream.pressure*pressureRatio;
rho2=freestream.density*densityRatio;
T2=freestream.temperature*temperatureRatio;
U2=freestream.speed/densityRatio;
p02=p2*(1+0.5*(gamma-1)*M2^2)^(gamma/(gamma-1));
T02=T2*(1+0.5*(gamma-1)*M2^2);
cp=gamma*freestream.R/(gamma-1);

%% 3. 保留守恒残差、方法来源与不能推出的激波几何
mass1=freestream.density*freestream.speed;
mass2=rho2*U2;
momentum1=freestream.pressure+freestream.density*freestream.speed^2;
momentum2=p2+rho2*U2^2;
energy1=cp*freestream.temperature+0.5*freestream.speed^2;
energy2=cp*T2+0.5*U2^2;
state=struct;
state.available=true;
state.method='完全气体法向激波 Rankine-Hugoniot 跃迁';
state.scope='小攻角球形钝鼻轴向驻点前的局部中心线状态；不求脱体距离或激波面形状';
state.upstreamMach=M1;
state.postShockMach=M2;
state.postShockPressure=p2;
state.postShockTemperature=T2;
state.postShockDensity=rho2;
state.postShockSpeed=U2;
state.postShockTotalPressure=p02;
state.postShockTotalTemperature=T02;
state.stagnationPressure=p02;
state.stagnationTemperature=T02;
state.pressureRatio=pressureRatio;
state.densityRatio=densityRatio;
state.temperatureRatio=temperatureRatio;
state.totalPressureRatio=p02/freestream.totalPressure;
state.entropyIncrease=freestream.R*log(freestream.totalPressure/p02);
state.shockStandOffDistance=NaN;
state.shockSurfacePosition=[NaN NaN NaN];
state.massResidual=(mass2-mass1)/mass1;
state.momentumResidual=(momentum2-momentum1)/momentum1;
state.energyResidual=(energy2-energy1)/energy1;
state.status='局部法向激波后状态已计算；脱体激波几何仍未求解';
end
