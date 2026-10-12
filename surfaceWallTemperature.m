function temperature = surfaceWallTemperature(faceCount,cfg)
%SURFACEWALLTEMPERATURE 读取面元壁温场，统一标量与非均匀壁温接口。
% 模块作用：cfg.flight.Twall 保留为初始/驻点单点壁温；可选
% cfg.heating.wallTemperatureByFace 为当前面元壁温 [K]，必须与网格
% 一一对应。NaN 明确表示该面元壁温未知，不用平均值或初始值补齐。
% 输出始终为 N×1，避免行向量隐式广播产生 N×N 的假热流矩阵。
validateattributes(faceCount,{'numeric'},{'scalar','integer','positive'});
if ~isfield(cfg,'flight') || ~isfield(cfg.flight,'Twall')
    error('RocketAeroHeat:WallField:Initial','缺少 cfg.flight.Twall。');
end
validateattributes(cfg.flight.Twall,{'numeric'}, ...
    {'real','finite','scalar','positive'});
temperature=repmat(double(cfg.flight.Twall),faceCount,1);
if isfield(cfg,'heating') && isfield(cfg.heating,'wallTemperatureByFace')
    value=cfg.heating.wallTemperatureByFace;
    validateattributes(value,{'numeric'},{'real','vector','numel',faceCount});
    value=double(value(:));
    if any(isinf(value) | value<=0)
        error('RocketAeroHeat:WallField:Temperature', ...
            '面元壁温须为正数；未知壁温用 NaN，不允许 Inf 或非正值。');
    end
    temperature=value;
end
end
