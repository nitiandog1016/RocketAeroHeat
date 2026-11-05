function value = idealAirThermodynamicState(mode,first,second,eddyViscosity)
%IDEALAIRTHERMODYNAMICSTATE 固定组分理想空气的显式变比热状态入口。
% 模块用途。采用Lemmon2000理想Helmholtz项，固定CoolProp7.2.0 OTH参考。
% 支持rhoT、PT、rhoe、Prho、hs、Ps，输入依次对应状态对字母，全部使用SI。
% e为热内能，调用者必须从总能量中先扣除动能及SST湍动能k。
% 数值温区60—2000K，超域或非法状态整次报错。输入同尺寸或一端标量。
% 输出T/rho/p/e/h/s/cp/cv/gamma/sound、指定自变量下偏导和工程输运。
% 第四参数为给定湍黏度mu_t，默认0，导数在mu_t冻结条件下计算。
% 输运采用Sutherland及Pr=.72、Pr_t=.9，其适用性须随工况另外验证。

%% 1. 读取固定系数并核对输入，防止无意的矩阵隐式广播
persistent data
if isempty(data)
    path = fullfile(fileparts(mfilename('fullpath')),'third_party','ideal_air','air_ideal.json');
    data = jsondecode(fileread(path));
    assert(strcmp(data.model,'Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'), ...
        'RocketAeroHeat:IdealAirData','热物性数据版本不匹配。');
end
mode = validatestring(mode,{'rhoT','PT','rhoe','Prho','hs','Ps'},mfilename,'mode');
validateattributes(first,{'numeric'},{'real','finite','nonempty'},mfilename,'first');
validateattributes(second,{'numeric'},{'real','finite','nonempty'},mfilename,'second');
[first,second] = matchingArrays(double(first),double(second));
if nargin < 4, eddyViscosity = 0; end
validateattributes(eddyViscosity,{'numeric'},{'real','finite','nonnegative','nonempty'});
[eddyViscosity,~] = matchingArrays(double(eddyViscosity),first);
assert(isequal(size(eddyViscosity),size(first)), ...
    'RocketAeroHeat:IdealAirShape','湍黏度应为标量或与状态同尺寸。');
R = data.gasConstant;

%% 2. 根据状态对恢复温度和密度，反解采用单调有界迭代
switch mode
    case 'rhoT'
        rho = first; temperature = second;
    case 'PT'
        temperature = second; rho = first./(R*temperature);
    case 'rhoe'
        rho = first; temperature = inverseTemperature(second,'e',ones(size(first))*101325,data);
    case 'Prho'
        assert(all(first(:)>0 & second(:)>0),'RocketAeroHeat:IdealAirState','压力和密度必须为正。');
        rho = second; temperature = first./(R*rho);
        % 仅对由p/rho运算产生的端点舍入作至多8eps量级恢复。
        tolerance = 8*eps*data.temperatureRangeK(2);
        assert(all(temperature(:)>=data.temperatureRangeK(1)-tolerance & ...
            temperature(:)<=data.temperatureRangeK(2)+tolerance), ...
            'RocketAeroHeat:IdealAirRange','反算温度超出60—2000K。');
        temperature = min(data.temperatureRangeK(2),max(data.temperatureRangeK(1),temperature));
    case 'hs'
        temperature = inverseTemperature(first,'h',ones(size(first))*101325,data);
        reference = assemble(ones(size(first))*data.reducingDensity,temperature,data);
        rho = data.reducingDensity*exp((reference.s-second)/R);
    case 'Ps'
        assert(all(first(:)>0),'RocketAeroHeat:IdealAirState','压力必须为正。');
        temperature = inverseTemperature(second,'s',first,data);
        rho = first./(R*temperature);
end
value = assemble(rho,temperature,data);

%% 3. 独立标识工程输运及温度偏导，导热系数包含比热变化项
value.mu = 1.716e-5*(temperature/273.15).^1.5*383.55./(temperature+110.4);
value.dmu_dt = value.mu.*(1.5./temperature-1./(temperature+110.4));
factor = value.mu/.72+eddyViscosity/.9;
value.lambda = value.cp.*factor;
value.dlambda_dt = value.dcp_dt.*factor+value.cp.*value.dmu_dt/.72;
assert(all(isfinite(value.lambda(:)) & isfinite(value.dlambda_dt(:))), ...
    'RocketAeroHeat:IdealAirState','输运状态超出数值范围。');

%% 4. 返回来源与适用域，供以后统一工作台和求解器登记方法
value.model = data.model;
value.reference = data.reference;
value.sourceDOI = data.sourceDOI;
value.sourceSHA256 = data.sourceSHA256;
value.gasConstant = R;
value.temperatureRangeK = data.temperatureRangeK(:).';
value.transportModel = 'Sutherland-273.15K-1.716e-5Pa.s-S110.4K-Pr0.72-Prt0.9';
value.energyDefinition = 'thermal internal energy; kinetic energy and SST k excluded';
value.status = 'thermodynamic-state-only';
value.productionQualified = false;
end

function [first,second] = matchingArrays(first,second)
% 只允许同尺寸或标量扩展，避免行列向量意外生成二维状态网格。
if isscalar(first)
    first = repmat(first,size(second));
elseif isscalar(second)
    second = repmat(second,size(first));
else
    assert(isequal(size(first),size(second)), ...
        'RocketAeroHeat:IdealAirShape','状态输入应同尺寸或一端标量。');
end
end

function value = assemble(rho,temperature,data)
% 组装理想状态和偏导，所有偏导的固定变量写在字段名称中。
assert(all(isfinite(rho(:)) & rho(:)>0),'RocketAeroHeat:IdealAirState','密度必须为有限正值。');
assert(all(isfinite(temperature(:)) & temperature(:)>=data.temperatureRangeK(1) & ...
    temperature(:)<=data.temperatureRangeK(2)), ...
    'RocketAeroHeat:IdealAirRange','温度超出60—2000K。');
[phi,first,second,third] = potential(temperature,data);
R = data.gasConstant;
tau = data.reducingTemperature./temperature;
value.T = temperature;
value.rho = rho;
value.p = rho.*R.*temperature;
value.e = R*data.reducingTemperature*first;
value.h = value.e+R*temperature;
value.s = R*(tau.*first-phi-log(rho/data.reducingDensity));
value.cv = -R*tau.^2.*second;
value.cp = value.cv+R;
value.gamma = value.cp./value.cv;
value.sound = sqrt(value.gamma.*R.*temperature);
value.dp_drho_e = R*temperature;
value.dp_de_rho = rho*R./value.cv;
value.dt_drho_e = zeros(size(rho));
value.dt_de_rho = 1./value.cv;
value.dh_drho_p = -value.cp.*temperature./rho;
value.dh_dp_rho = value.cp./(rho*R);
value.ds_drho_p = -value.cp./rho;
value.ds_dp_rho = value.cv./value.p;
value.dcp_dt = R./temperature.*(2*tau.^2.*second+tau.^3.*third);
fields = fieldnames(value);
for index = 1:numel(fields)
    field = value.(fields{index});
    assert(all(isfinite(field(:))),'RocketAeroHeat:IdealAirState','状态超出数值范围。');
end
assert(all(value.cv(:)>0 & value.p(:)>0),'RocketAeroHeat:IdealAirState','状态失去正性。');
end

function [phi,first,second,third] = potential(temperature,data)
% 对固定理想Helmholtz项解析求导，独立于原NASA7工程模型。
terms = data.alpha0;
power = terms{2}; logarithm = terms{3}; oscillator = terms{4};
generalized = terms{5}; offset = terms{6};
tau = data.reducingTemperature./temperature;
phi = offset.a1+offset.a2*tau;
first = ones(size(tau))*offset.a2;
second = zeros(size(tau)); third = second;
for index = 1:numel(power.n)
    exponent = power.t(index);
    term = power.n(index)*tau.^exponent;
    phi = phi+term;
    first = first+term*exponent./tau;
    second = second+term*exponent*(exponent-1)./tau.^2;
    third = third+term*exponent*(exponent-1)*(exponent-2)./tau.^3;
end
phi = phi+logarithm.a*log(tau);
first = first+logarithm.a./tau;
second = second-logarithm.a./tau.^2;
third = third+2*logarithm.a./tau.^3;
for index = 1:numel(oscillator.n)
    frequency = oscillator.t(index); weight = oscillator.n(index);
    decay = exp(-frequency*tau); complement = -expm1(-frequency*tau);
    phi = phi+weight*log(complement);
    first = first+weight*frequency*decay./complement;
    second = second-weight*frequency^2*decay./complement.^2;
    third = third+weight*frequency^3*decay.*(1+decay)./complement.^3;
end
weight = generalized.n; frequency = generalized.t;
decay = generalized.c*exp(-frequency*tau); denominator = 1+decay;
phi = phi+weight*(frequency*tau+log1p(decay));
first = first+weight*frequency./denominator;
second = second+weight*frequency^2*decay./denominator.^2;
third = third-weight*frequency^3*decay.*(1-decay)./denominator.^3;
end

function temperature = inverseTemperature(target,quantity,pressure,data)
% 每个状态保持独立温度括区间，Newton越界则二分，至多60次迭代。
lower = ones(size(target))*data.temperatureRangeK(1);
upper = ones(size(target))*data.temperatureRangeK(2);
low = assemble(pressure./(data.gasConstant*lower),lower,data);
high = assemble(pressure./(data.gasConstant*upper),upper,data);
minimum = low.(quantity); maximum = high.(quantity);
assert(all(target(:)>=minimum(:) & target(:)<=maximum(:)), ...
    'RocketAeroHeat:IdealAirRange','反解目标超出60—2000K状态范围。');
temperature = lower+(upper-lower).*(target-minimum)./(maximum-minimum);
temperature(target==minimum) = data.temperatureRangeK(1);
temperature(target==maximum) = data.temperatureRangeK(2);
for iteration = 1:60
    value = assemble(pressure./(data.gasConstant*temperature),temperature,data);
    residual = value.(quantity)-target;
    active = abs(residual)>2e-13*max(abs(target),1);
    if ~any(active(:)), return; end
    switch quantity
        case 'e', derivative = value.cv;
        case 'h', derivative = value.cp;
        case 's', derivative = value.cp./temperature;
    end
    raiseLower = active & residual<0;
    lower(raiseLower) = temperature(raiseLower);
    lowerUpper = active & residual>=0;
    upper(lowerUpper) = temperature(lowerUpper);
    candidate = temperature-residual./derivative;
    bisect = candidate<=lower | candidate>=upper;
    candidate(bisect) = (lower(bisect)+upper(bisect))/2;
    temperature(active) = candidate(active);
end
error('RocketAeroHeat:IdealAirInverse','温度反解在60次迭代内未收敛。');
end
