function total=surfaceFlowTotalStateFromAnchor(anchor,gas)
%SURFACEFLOWTOTALSTATEFROMANCHOR 从实际物理锚点显式转换热力总状态约定。
% 模块一：输入p/T/rho/U保持原值；不能把旧总压直接套到新熵积分常数。
% 模块二：同一cp积分焓求T0，再由同一标准熵求p0，只转换一次。
% 模块三：保存原总状态及差量；转换不代表新求物理外流或新增覆盖。
required={'pressure','temperature','density','speed','method'};
assert(isstruct(anchor)&&isscalar(anchor)&&all(isfield(anchor,required)), ...
    'SurfaceFlow:Anchor','需要实际锚点p/T/rho/U及来源方法。');
if isfield(anchor,'available'),assert(isequal(anchor.available,true),'SurfaceFlow:Anchor','锚点求解未成功。');end
for name={'pressure','temperature','density','speed'}
    validateattributes(anchor.(name{1}),{'double'},{'real','finite','scalar','positive'});
end
state=surfaceFlowThermodynamicState(anchor.temperature,anchor.pressure,gas);
assert(abs(anchor.pressure-anchor.density*gas.R*anchor.temperature)<= ...
    1e-9*max(anchor.pressure,anchor.density*gas.R*anchor.temperature), ...
    'SurfaceFlow:EOS','锚点状态方程不一致。');
target=state.enthalpy+anchor.speed^2/2;
upper=surfaceFlowThermodynamicState(gas.temperatureRangeK(2),anchor.pressure,gas);
assert(target<=upper.enthalpy,'SurfaceFlow:Temperature','锚点总焓超出声明温域。');
lo=anchor.temperature;hi=gas.temperatureRangeK(2);
for k=1:47
    mid=(lo+hi)/2;value=surfaceFlowThermodynamicState(mid,anchor.pressure,gas);
    if value.enthalpy<target,lo=mid;else,hi=mid;end
end
T0=(lo+hi)/2;atT0=surfaceFlowThermodynamicState(T0,anchor.pressure,gas);
p0=anchor.pressure*exp((atT0.entropy-state.entropy)/gas.R);
assert(isfinite(p0)&&p0>anchor.pressure,'SurfaceFlow:Anchor','转换后的总压必须有限且高于静压。');
total=struct('available',true,'totalPressure',p0,'totalTemperature',T0, ...
    'thermalPropertyMode',state.thermalPropertyMode,'thermalPropertyConvention',gas.thermalPropertyConvention, ...
    'method',[char(anchor.method) '；从原物理p/T/rho/U转换为连续cp积分总状态'], ...
    'sourceAnchor',struct('pressure',anchor.pressure,'temperature',anchor.temperature, ...
    'density',anchor.density,'speed',anchor.speed,'method',anchor.method), ...
    'sourcePhysicalStateChanged',false,'totalPressureResetAlongRay',false, ...
    'newExternalPressureOrHeatingSolved',false);
if isfield(anchor,'totalPressure')&&isfield(anchor,'totalTemperature')
    total.originalTotalPressure=anchor.totalPressure;total.originalTotalTemperature=anchor.totalTemperature;
    total.totalPressureConventionDifference=p0-anchor.totalPressure;
    total.totalTemperatureConventionDifference=T0-anchor.totalTemperature;
end
end
