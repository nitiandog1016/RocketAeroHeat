function boundary=wallBackBoundary(material)
%WALLBACKBOUNDARY 验证壁面内侧边界；省略 backBoundary 时保持绝热。
% 模块一：识别 type='adiabatic' 或 'convectionRadiation'。
% 模块二：校验用户给定的恒定内侧 h、流体温度、有效辐射率与辐射温度。
% 有效辐射率包含用户已知的交换条件；这里不求内腔流动、视角因子或
% 封闭腔体的辐射网络。不能从 STL 自动猜测内环境参数。
boundary=struct('type','adiabatic','active',false,'heatTransferCoefficient',0, ...
    'fluidTemperature',0,'emissivity',0,'radiationTemperature',0, ...
    'method','背面绝热：内侧热流为零', ...
    'scope','平壁等面积近似；背面无能量交换');
if ~isfield(material,'backBoundary') || isempty(material.backBoundary), return; end
input=material.backBoundary;
if ~isstruct(input) || ~isscalar(input) || ~isfield(input,'type') || ...
        ~(ischar(input.type) || (isstring(input.type) && isscalar(input.type)))
    error('RocketAeroHeat:BackBoundary:Type','backBoundary 需要标量结构与明确的 type。');
end
kind=char(input.type);
if strcmp(kind,'adiabatic'), return; end
if ~strcmp(kind,'convectionRadiation')
    error('RocketAeroHeat:BackBoundary:Type','type 须为 adiabatic 或 convectionRadiation。');
end
for key={'heatTransferCoefficient','fluidTemperature','emissivity','radiationTemperature'}
    if ~isfield(input,key{1})
        error('RocketAeroHeat:BackBoundary:Missing','内侧边界缺少 %s。',key{1});
    end
    validateattributes(input.(key{1}),{'numeric'},{'real','finite','scalar'});
    boundary.(key{1})=double(input.(key{1}));
end
if boundary.heatTransferCoefficient<0 || boundary.fluidTemperature<=0 || ...
        boundary.emissivity<0 || boundary.emissivity>1 || boundary.radiationTemperature<=0
    error('RocketAeroHeat:BackBoundary:Range','内侧 h>=0、两个环境温度>0，有效辐射率在 0~1。');
end
boundary.type=kind; boundary.active=true;
boundary.method='背面 Newton 对流 + 灰体辐射（指定恒定内环境）';
boundary.scope=['内侧 h、流体温度和辐射环境为用户指定的均匀常数；' ...
    '内外面积相同的平壁近似，不求内腔流动或封闭腔体辐射网络'];
end
