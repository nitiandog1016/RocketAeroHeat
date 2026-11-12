function state=phaseAwareAirReferenceState(p,Te,Ue,Tw,recoveryFactor,transportModel)
%PHASEAWAREAIRREFERENCESTATE 同一标准空气的恢复焓、参考焓及输运接口。
% 模块用途：接收显式外缘状态和有来源的恢复因子r，计算haw=he+rUe²/2，
% h*=0.5(he+hw)+0.22(haw-he)。只提供相容热接口，不自动赋予外缘或热流资格。
% 输入为同尺寸或标量，Ue为速度大小。r的层流/湍流来源由调用区域登记。
% 第六项显式选择解析或线性表输运；省略时保留原解析实现。
if nargin<6, transportModel='lemmonAnalytic'; end

%% 1. 确认输入形状与恢复因子，拒绝隐式二维广播
values={p,Te,Ue,Tw,recoveryFactor}; sizes=cellfun(@size,values,'UniformOutput',false);
nonScalar=find(cellfun(@numel,values)>1); shape=[1 1];
if ~isempty(nonScalar), shape=sizes{nonScalar(1)}; end
for k=1:5
    validateattributes(values{k},{'numeric'},{'real','finite','nonempty'});
    if isscalar(values{k}), values{k}=repmat(double(values{k}),shape);
    else, assert(isequal(size(values{k}),shape),'RAH:PhaseReference:Shape','输入尺寸不一致。'); end
end
[p,Te,Ue,Tw,r]=deal(values{:});
assert(all(Ue(:)>=0&r(:)>=0&r(:)<=1),'RAH:PhaseReference:Recovery','速度或恢复因子无效。');

%% 2. 同模型能量运算，不把NASA7/OTH参考值静默传入
edge=phaseAwareAirState('PT',p,Te,transportModel); wall=phaseAwareAirState('PT',p,Tw,transportModel);
haw=edge.enthalpyJkg+r.*Ue.^2/2; recovered=phaseAwareAirState('Ph',p,haw,transportModel);
href=.5*(edge.enthalpyJkg+wall.enthalpyJkg)+.22*(haw-edge.enthalpyJkg);
reference=phaseAwareAirState('Ph',p,href,transportModel);

%% 3. 返回区分身份的温度和壁温导数；这里没有求固体壁温或壁热流
state=struct('model',edge.model,'domain',edge.domain,'sourceSHA256',edge.sourceSHA256, ...
    'method','Riley NASA-TM-107838 Appendix H.6/H.7 + Lemmon standard-air Cp', ...
    'givenWallTemperatureK',Tw,'edgeTemperatureInputK',Te, ...
    'recoveryTemperatureK',recovered.temperatureK,'referenceTemperatureK',reference.temperatureK, ...
    'recoveryEnthalpyJkg',haw,'referenceEnthalpyJkg',href, ...
    'thermalDrivingEnthalpyJkg',haw-wall.enthalpyJkg, ...
    'referenceTemperatureDerivativeByWallTemperature',.5*wall.cpJkgK./reference.cpJkgK, ...
    'dynamicViscosityPas',reference.dynamicViscosityPas,'conductivityWmK',reference.conductivityWmK, ...
    'Prandtl',reference.Prandtl,'recoveryFactorInput',r, ...
    'referenceDensityKgm3',reference.densityKgm3,'transportModel',reference.transportModel, ...
    'wallHeatFluxWm2',[],'solidWallTemperatureK',[],'physicalQualification',false,'engineeringHeatReady',false);
end
