function freestream = calcFreestream(flight)
%% ========================================================================
% calcFreestream.m
%
% Rocket AeroHeat
%
% 三维自由来流状态计算模块
%
% ========================================================================
%
% 本文件是纯计算模块。
%
% 不创建 GUI。
% 不绘图。
% 不修改 mesh。
%
%
% 输入：
%
%       flight
%
%
% 必需字段：
%
%       flight.altitude      [m]
%
%       flight.Mach          [-]
%
%       flight.alpha         [deg]
%
%       flight.beta          [deg]
%
%
% 可选字段：
%
%       flight.atmosphere
%
% 如果不存在，则自动调用：
%
%       standardAtmosphere()
%
%
% ========================================================================
%
% 坐标约定
%
%       +X = 箭体前向 / 鼻头方向
%
%       +Y = 侧向
%
%       +Z = 向上
%
%
% 零攻角、零侧滑时：
%
%       Vehicle direction:
%
%           +X
%
%
%       Air flow relative to vehicle:
%
%           -X
%
%
% ========================================================================
%
% 攻角 alpha 的约定
%
% 正 alpha：
%
%       表示箭体鼻部相对于飞行轨迹抬头。
%
%
% 因此在箭体坐标系中：
%
%       vehicle velocity
%
% 会向 -Z 偏转。
%
%
% 相对空气来流：
%
%       flow direction
%
% 则向 +Z 偏转。
%
%
% ========================================================================
%
% 侧滑角 beta 的约定
%
% 正 beta：
%
%       vehicle velocity 向 +Y 偏转。
%
%
% 因此相对来流：
%
%       flow direction 向 -Y 偏转。
%
%
% ========================================================================
%
% 最终方向：
%
% Vehicle unit direction：
%
%   [ cos(alpha)cos(beta),
%     sin(beta),
%    -sin(alpha)cos(beta) ]
%
%
% Flow unit direction：
%
%       -Vehicle unit direction
%
%
% ========================================================================
%
% 后续迎风判据统一使用：
%
%       mu = -n dot flowDirection
%
%
% 当：
%
%       mu > 0
%
% 为迎风。
%
%
%       mu <= 0
%
% 为背风。
%
%
%% ========================================================================



%% ========================================================================
% 1. 输入检查
% ========================================================================

if nargin < 1 || ~isstruct(flight) || ~isscalar(flight)

    error( ...
        'calcFreestream: 必须输入 flight 结构体。');

end


requiredFields = { ...
    'altitude', ...
    'Mach', ...
    'alpha', ...
    'beta'};


for i = 1:numel(requiredFields)

    fieldName = ...
        requiredFields{i};


    if ~isfield(flight,fieldName)

        error( ...
            'calcFreestream: 缺少 flight.%s。', ...
            fieldName);

    end

    validateattributes(flight.(fieldName), {'numeric'}, ...
        {'real','finite','scalar'}, mfilename, ['flight.' fieldName]);

end



%% ========================================================================
% 2. 检查输入是否为有限标量
% ========================================================================

inputValues = [ ...
    flight.altitude, ...
    flight.Mach, ...
    flight.alpha, ...
    flight.beta];


if any(~isfinite(inputValues)) || ...
        any(~isreal(inputValues))

    error( ...
        'calcFreestream: flight 中存在非法数值。');

end


if flight.Mach < 0

    error( ...
        'calcFreestream: Mach 数不能小于 0。');

end



%% ========================================================================
% 3. 获取标准大气状态
% ========================================================================

if isfield(flight,'atmosphere') && ...
        ~isempty(flight.atmosphere)

    atmosphere = ...
        flight.atmosphere;

    % Cached standard-atmosphere data must match the requested altitude.
    if isstruct(atmosphere) && isscalar(atmosphere) && ...
            isfield(atmosphere,'geometricAltitude') && ...
            ~isequal(double(atmosphere.geometricAltitude), double(flight.altitude))
        atmosphere = standardAtmosphere(flight.altitude);
    end

else

    atmosphere = ...
        standardAtmosphere( ...
            flight.altitude);

end



%% ========================================================================
% 4. 基本物理量
% ========================================================================

validateattributes(flight.altitude, {'numeric'}, {'nonnegative'});
if ~isstruct(atmosphere) || ~isscalar(atmosphere)
    error('calcFreestream:InvalidAtmosphere', 'atmosphere must be a scalar struct.');
end
physicalFields = {'gamma','R','temperature','pressure','density', ...
    'speedOfSound','dynamicViscosity'};
for i = 1:numel(physicalFields)
    name = physicalFields{i};
    if ~isfield(atmosphere,name)
        error('calcFreestream:InvalidAtmosphere', 'Missing atmosphere.%s.',name);
    end
    validateattributes(atmosphere.(name), {'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['atmosphere.' name]);
end
validateattributes(atmosphere.gamma, {'numeric'}, {'>',1});

Mach = ...
    flight.Mach;


alphaDeg = ...
    flight.alpha;


betaDeg = ...
    flight.beta;


alpha = ...
    deg2rad(alphaDeg);


beta = ...
    deg2rad(betaDeg);


gamma = ...
    atmosphere.gamma;


R = ...
    atmosphere.R;


Tinf = ...
    atmosphere.temperature;


pinf = ...
    atmosphere.pressure;


rhoinf = ...
    atmosphere.density;


ainf = ...
    atmosphere.speedOfSound;


muinf = ...
    atmosphere.dynamicViscosity;



%% ========================================================================
% 5. 自由来流速度大小
%
%       V∞ = M∞ a∞
%% ========================================================================

speed = ...
    Mach ...
    * ...
    ainf;



%% ========================================================================
% 6. 飞行方向单位向量
% ========================================================================
%
% 坐标：
%
%       X 前
%       Y 侧
%       Z 上
%
%
% alpha > 0：
%
%       vehicle velocity 朝 -Z
%
%
% beta > 0：
%
%       vehicle velocity 朝 +Y
%

ca = cos(alpha);

sa = sin(alpha);

cb = cos(beta);

sb = sin(beta);


unitFlightDirection = [ ...
    ca*cb, ...
    sb, ...
   -sa*cb];


%% ------------------------------------------------------------------------
% 数值保险：重新单位化
%% ------------------------------------------------------------------------

unitFlightDirection = ...
    unitFlightDirection ...
    / ...
    max( ...
        norm(unitFlightDirection), ...
        eps);



%% ========================================================================
% 7. 空气相对箭体的流动方向
% ========================================================================
%
% 空气相对箭体运动方向：
%
%       与箭体相对空气速度相反。
%

unitFlowDirection = ...
    -unitFlightDirection;



%% ========================================================================
% 8. 飞行速度矢量
% ========================================================================

flightVelocityVector = ...
    speed ...
    * ...
    unitFlightDirection;



%% ========================================================================
% 9. 空气相对箭体速度矢量
% ========================================================================

flowVelocityVector = ...
    speed ...
    * ...
    unitFlowDirection;



%% ========================================================================
% 10. 动压
%
%             1
%       q∞ = --- rho∞ V∞²
%             2
%% ========================================================================

dynamicPressure = ...
    0.5 ...
    * ...
    rhoinf ...
    * ...
    speed^2;



%% ========================================================================
% 11. 单位长度 Reynolds 数
%
%             rho∞ V∞
% Re / L = -------------
%                mu∞
%
% 单位：
%
%       1/m
%% ========================================================================

if muinf > 0

    reynoldsPerMeter = ...
        rhoinf ...
        * ...
        speed ...
        / ...
        muinf;

else

    reynoldsPerMeter = ...
        NaN;

end



%% ========================================================================
% 12. 理想气体等熵总温
%
%                         gamma - 1
% T0 = T∞ [ 1 + ---------------- M² ]
%                             2
%
% 注意：
%
% 这是来流总状态。
%
% 不是钝头激波之后的驻点压力。
%% ========================================================================

totalTemperature = ...
    Tinf ...
    * ...
    ( ...
    1 ...
    + ...
    0.5 ...
    * ...
    (gamma-1) ...
    * ...
    Mach^2);



%% ========================================================================
% 13. 理想气体等熵总压
%% ========================================================================

totalPressure = ...
    pinf ...
    * ...
    ( ...
    1 ...
    + ...
    0.5 ...
    * ...
    (gamma-1) ...
    * ...
    Mach^2 ...
    ) ...
    ^ ...
    ( ...
    gamma ...
    / ...
    (gamma-1));



%% ========================================================================
% 14. 单位动能
%% ========================================================================

kineticEnergyPerMass = ...
    0.5 ...
    * ...
    speed^2;



%% ========================================================================
% 15. 输出结构体
% ========================================================================

freestream = struct;



%% ------------------------------------------------------------------------
% 飞行条件
%% ------------------------------------------------------------------------

freestream.altitude = ...
    flight.altitude;


freestream.Mach = ...
    Mach;


freestream.alpha = ...
    alphaDeg;


freestream.beta = ...
    betaDeg;


freestream.alphaRad = ...
    alpha;


freestream.betaRad = ...
    beta;



%% ------------------------------------------------------------------------
% 自由来流热力状态
%% ------------------------------------------------------------------------

freestream.temperature = ...
    Tinf;


freestream.pressure = ...
    pinf;


freestream.density = ...
    rhoinf;


freestream.speedOfSound = ...
    ainf;


freestream.dynamicViscosity = ...
    muinf;


freestream.gamma = ...
    gamma;


freestream.R = ...
    R;



%% ------------------------------------------------------------------------
% 速度
%% ------------------------------------------------------------------------

freestream.speed = ...
    speed;


freestream.Vinf = ...
    speed;


freestream.unitFlightDirection = ...
    unitFlightDirection;


freestream.unitFlowDirection = ...
    unitFlowDirection;


freestream.flightVelocityVector = ...
    flightVelocityVector;


freestream.flowVelocityVector = ...
    flowVelocityVector;



%% ------------------------------------------------------------------------
% 动压 / Reynolds 数
%% ------------------------------------------------------------------------

freestream.dynamicPressure = ...
    dynamicPressure;


freestream.qinf = ...
    dynamicPressure;


freestream.reynoldsPerMeter = ...
    reynoldsPerMeter;



%% ------------------------------------------------------------------------
% 总状态
%% ------------------------------------------------------------------------

freestream.totalTemperature = ...
    totalTemperature;


freestream.totalPressure = ...
    totalPressure;


freestream.kineticEnergyPerMass = ...
    kineticEnergyPerMass;



%% ------------------------------------------------------------------------
% 标准大气
%% ------------------------------------------------------------------------

freestream.atmosphere = ...
    atmosphere;



%% ------------------------------------------------------------------------
% 坐标定义
%% ------------------------------------------------------------------------

freestream.coordinateConvention.bodyXAxis = ...
    [1 0 0];


freestream.coordinateConvention.bodyYAxis = ...
    [0 1 0];


freestream.coordinateConvention.bodyZAxis = ...
    [0 0 1];


freestream.coordinateConvention.description = ...
    ['+X forward/nose, +Y lateral, +Z upward. ', ...
     'Positive alpha gives flow toward +Z. ', ...
     'Positive beta gives flow toward -Y.'];



%% ------------------------------------------------------------------------
% 后续面元判据
%% ------------------------------------------------------------------------

freestream.panelConvention.windwardMetric = ...
    'mu = -dot(n, unitFlowDirection)';


freestream.panelConvention.windwardCondition = ...
    'mu > 0';


freestream.panelConvention.leewardCondition = ...
    'mu <= 0';

end
