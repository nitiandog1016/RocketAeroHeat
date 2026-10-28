function result = solveAxisymmetricBaseCorrelation(fs,shoulder,geometry,Twall,cfg)
%SOLVEAXISYMMETRICBASECORRELATION 轴对称无喷流底压与代表底热关联。
% 模块用途：以肩部已经求得的外流状态及边界层类别，复算
% Lamb / Oberkampf, SAND93-0280 (1993) 的轴对称底部关联。
% 输入全部采用 SI 单位；主轴 X，鼻尖 +X，零攻角来流 -X。
%
% fs：pressure, temperature, density, speed, Mach, gamma, R, alpha, beta。
% shoulder：pressure, temperature, density, speed, Mach,
%           dynamicViscosity, distance，以及逐输入 sources.* 字符说明。
% geometry：family='cone'/'cylinder', boundaryLayerRegime='laminar'/'turbulent',
%           baseRadius, noseBluntnessRatio, hasRearSupport, hasGasInjection,
%           source；圆锥还需要 halfAngleDeg。
% Twall：给定壁温 [K]，不是本模块求解的壁温。
%
% 圆锥层流：报告 Eq.(8b)、Eq.(17)；圆锥湍流：Eq.(11)、Eq.(18)。
% 圆柱湍流只有底压 Eq.(10)，不套用圆锥的底热公式。
% St_b=q_b/[rho_inf U_inf (H_inf-h_w)]，使用自由流总焓差。
% 特别注意 Eq.(18) 的整个压力组合先开平方再取 1.85 次方。
%
% 输出代表底热量；报告未提供可核验的径向面分布/面积平均定义，
% 因此不能把该标量涂满整个底面或计入已求解热流面积。底部外缘
% 静温、回流速度与求解壁温保持 NaN。不同输入来源单独保留。

%% 1. 输入结构及逐字段来源检查，建立真实缺失值输出
if nargin<5, cfg=struct; end
fields={'pressure','temperature','density','speed','Mach','gamma','R','alpha','beta'};
for k=1:numel(fields)
    key=fields{k};
    if ~isfield(fs,key)
        error('RocketAeroHeat:BaseCorrelation:Input','缺少 fs.%s。',key);
    end
    validateattributes(fs.(key),{'numeric'},{'real','finite','scalar'});
end
fields={'pressure','temperature','density','speed','Mach','dynamicViscosity','distance'};
for k=1:numel(fields)
    key=fields{k};
    if ~isfield(shoulder,key) || ~isfield(shoulder,'sources') || ...
            ~isfield(shoulder.sources,key) || isempty(shoulder.sources.(key))
        error('RocketAeroHeat:BaseCorrelation:ShoulderSource', ...
            '肩部输入 %s 及其 sources.%s 必须同时提供。',key,key);
    end
    validateattributes(shoulder.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'});
end
for key={'family','boundaryLayerRegime','baseRadius','noseBluntnessRatio', ...
        'hasRearSupport','hasGasInjection','flatCircularBase','isolatedAxisymmetricBody','source'}
    if ~isfield(geometry,key{1})
        error('RocketAeroHeat:BaseCorrelation:Geometry','缺少 geometry.%s。',key{1});
    end
end
validateattributes(geometry.baseRadius,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(geometry.noseBluntnessRatio,{'numeric'}, ...
    {'real','finite','scalar','>=',0});
validateattributes(geometry.hasRearSupport,{'logical','numeric'},{'binary','scalar'});
validateattributes(geometry.hasGasInjection,{'logical','numeric'},{'binary','scalar'});
validateattributes(geometry.flatCircularBase,{'logical','numeric'},{'binary','scalar'});
validateattributes(geometry.isolatedAxisymmetricBody,{'logical','numeric'},{'binary','scalar'});
validateattributes(Twall,{'numeric'},{'real','scalar'});
if isempty(geometry.source)
    error('RocketAeroHeat:BaseCorrelation:GeometrySource','几何来源不能为空。');
end
result=initialResult(fs,shoulder,geometry,Twall);
cp=fs.gamma*fs.R/(fs.gamma-1);
T0=fs.temperature+fs.speed^2/(2*cp);

%% 2. 严格限定已核验的气体、几何及流态范围
if fs.Mach<=1 || fs.gamma<=1 || fs.R<=0 || ...
        min([fs.pressure,fs.temperature,fs.density,fs.speed])<=0
    result.status='自由流物理状态无效'; return;
end
if abs(fs.alpha)>1e-8 || abs(fs.beta)>1e-8
    result.status='零攻角关联不适用于当前攻角或侧滑角'; return;
end
if abs(fs.gamma-1.4)>1e-8 || abs(fs.R/287.05287-1)>0.005
    result.status='本实现限 gamma=1.4 的定比热空气；未扩展真实气体'; return;
end
if geometry.hasRearSupport || geometry.hasGasInjection
    result.status='尾撑或底部喷流改变回流区，不适用无干扰底部关联'; return;
end
if ~geometry.flatCircularBase || ~geometry.isolatedAxisymmetricBody
    result.status='本关联只适用独立轴对称物体的平圆底面'; return;
end
if geometry.noseBluntnessRatio>0.3
    result.status='超出本实现的细长小钝度几何范围 rn/rb<=0.3'; return;
end
relativeResidual=max([abs(fs.density*fs.R*fs.temperature/fs.pressure-1), ...
    abs(fs.speed/(fs.Mach*sqrt(fs.gamma*fs.R*fs.temperature))-1), ...
    abs(shoulder.density*fs.R*shoulder.temperature/shoulder.pressure-1), ...
    abs(shoulder.speed/(shoulder.Mach*sqrt(fs.gamma*fs.R*shoulder.temperature))-1), ...
    abs((cp*shoulder.temperature+shoulder.speed^2/2)/(cp*T0)-1)]);
result.inputStateRelativeResidual=relativeResidual;
if relativeResidual>1e-5
    result.status='自由流/肩部理想气体及共同总焓不一致'; return;
end
if fs.temperature<80 || shoulder.temperature<80 || T0>2500
    result.status='超出本实现输运/无化学空气筛选范围 80<=T<=2500 K'; return;
end
isCone=strcmp(geometry.family,'cone');
isCylinder=strcmp(geometry.family,'cylinder');
isLaminar=strcmp(geometry.boundaryLayerRegime,'laminar');
isTurbulent=strcmp(geometry.boundaryLayerRegime,'turbulent');
if ~(isCone||isCylinder) || ~(isLaminar||isTurbulent)
    result.status='仅支持明确的圆锥/圆柱与层流/湍流类别'; return;
end
if isCone
    if ~isfield(geometry,'halfAngleDeg')
        error('RocketAeroHeat:BaseCorrelation:ConeAngle','圆锥需要 halfAngleDeg。');
    end
    angle=geometry.halfAngleDeg;
    validateattributes(angle,{'numeric'},{'real','finite','scalar','positive'});
    if (isLaminar && (angle<9 || angle>10)) || ...
            (isTurbulent && (angle<6 || angle>12))
        result.status=['本实现底热证据限层流9~10度、湍流6~12度圆锥；' ...
            '未采用报告对其他半锥角的临时外推']; return;
    end
    if fs.Mach<2 || fs.Mach>20 || shoulder.Mach<=1
        result.status='超出支持的超声速细长锥工况范围'; return;
    end
else
    if isLaminar
        result.status='本独立模块尚未建立圆柱层流底热或底压闭合'; return;
    end
    if shoulder.Mach<1.2 || shoulder.Mach>15
        result.status='圆柱底压关联要求肩部 1.2<=Mach<=15'; return;
    end
end
Re=shoulder.density*shoulder.speed*shoulder.distance/shoulder.dynamicViscosity;
result.shoulderReynolds=Re;
if isLaminar && Re<1e5
    result.status='报告指出层流 Re_e,s<1e5 时关联失效，不更改指数补值'; return;
end
knLimit=0.01;
if isfield(cfg,'maxContinuumKnudsen'), knLimit=cfg.maxContinuumKnudsen; end
validateattributes(knLimit,{'numeric'},{'real','finite','scalar','>',0,'<=',0.01});
muInf=1.458e-6*fs.temperature^1.5/(fs.temperature+110.4);
lambdaInf=airMeanFreePath(fs.pressure,fs.temperature,muInf,fs.R);
lambdaEdge=airMeanFreePath(shoulder.pressure,shoulder.temperature, ...
    shoulder.dynamicViscosity,fs.R);

%% 3. 计算底压；保留标准表面 Cp 符号，避免报告底阻系数符号混用
Me=shoulder.Mach; pe=shoulder.pressure;
if isCone && isLaminar
    C4=(1+2/(fs.gamma-1)/Me^2)^(-2);
    parameter=Me^2/sqrt(Re);
    ratio=3.05*parameter^0.6/C4;
    result.pressureMethod='SAND93-0280 Eq.(8b)：轴对称层流底压关联';
elseif isCone
    J=1.7/log(21/geometry.halfAngleDeg);
    ratio=(pe/fs.pressure)^2*(0.025+0.906/(1+(fs.gamma-1)/2*Me^2))^J;
    result.correlationExponentJ=J;
    result.pressureMethod='SAND93-0280 Eq.(11)：细长圆锥湍流底压关联';
else
    ratio=(pe/fs.pressure)*(0.05+0.967/(1+(fs.gamma-1)/2*Me^2));
    result.pressureMethod='SAND93-0280 Eq.(10)：圆柱湍流底压关联，以肩部为参考';
end
pb=fs.pressure*ratio;
% 背风底部应处于肩部静压以下；异常不做裁剪，也不伪装成附着解。
if ~isfinite(pb) || pb<=0 || pb>=pe
    result.status='关联给出的底压未形成低于肩部压力的回流区'; return;
end
% 低速底回流 Tb≈T0 的近似仅用于保守 Kn 筛选，不能当成底部静温解。
muBase=1.458e-6*T0^1.5/(T0+110.4);
lambdaBase=airMeanFreePath(pb,T0,muBase,fs.R);
result.continuumKnudsen=max([lambdaInf,lambdaEdge,lambdaBase])/geometry.baseRadius;
result.continuumMethod='自由流/肩部/底回流T0筛选近似的最大平均自由程除以底半径';
result.baseGasTemperatureUsedOnlyForKnScreen=T0;
if ~isfinite(result.continuumKnudsen) || result.continuumKnudsen>=knLimit
    result.status='底回流区 Kn 超出当前无滑移连续流关联范围'; return;
end
result.pressureAvailable=true;
result.basePressure=pb;
result.basePressureRatio=ratio;
result.basePressureCoefficient=(pb-fs.pressure)/(0.5*fs.density*fs.speed^2);
result.fieldMetadata.basePressure=fieldRecord('computed',result.pressureMethod, ...
    true,'代表底压；报告p11支持低回流速度下近似均匀静压，不含底唇细节');

%% 4. 代表底热关联：不生成没有依据的全底面径向分布
if isCylinder
    result.status='圆柱底压已计算；报告未提供相应圆柱底热关联';
    return;
end
if isLaminar
    normalizedParameter=ratio*C4;
    St=1e-4*normalizedParameter^1.8;
    result.heatingMethod='SAND93-0280 Eq.(17)：细长圆锥层流代表底Stanton关联';
else
    normalizedParameter=sqrt(ratio*(fs.pressure/pe)^2);
    St=8.3e-4*normalizedParameter^1.85;
    result.heatingMethod='SAND93-0280 Eq.(18)：6~12度圆锥湍流代表底Stanton关联';
end
result.baseStantonNumber=St;
result.heatingCorrelationParameter=normalizedParameter;
result.stantonAvailable=true;
result.fieldMetadata.baseStantonNumber=fieldRecord('computed', ...
    result.heatingMethod,true,result.heatingScope);
if ~isfinite(Twall) || Twall<80 || Twall>2500
    result.status='底压与代表St已计算；缺少适用范围内的给定壁温，热流保持NaN';
    return;
end
enthalpyDifference=cp*(T0-Twall);
q=St*fs.density*fs.speed*enthalpyDifference;
result.heatFluxAvailable=true;
result.heatFluxRepresentative=q;
result.heatTransferCoefficientRepresentative=St*fs.density*fs.speed*cp;
result.heatRecoveryReferenceTemperature=T0;
result.totalEnthalpyDifference=enthalpyDifference;
result.available=true;
result.status='底压与代表底热关联已计算；全底面的空间分布仍未建立';
result.fieldMetadata.heatFluxRepresentative=fieldRecord('computed', ...
    [result.heatingMethod '；q=St rho_inf U_inf (H_inf-h_w)'],true,result.heatingScope);
result.fieldMetadata.heatRecoveryReferenceTemperature=fieldRecord('computed', ...
    '报告p27底St定义的自由流总焓参考温度；不是底部静温/外缘温度',true, ...
    'gamma=1.4定比热空气的总能量参考');
end

function r=initialResult(fs,shoulder,geometry,Twall)
% 模块用途：初始化结果及来源，未求解量从始至终保持 NaN。
r=struct('available',false,'pressureAvailable',false,'stantonAvailable',false, ...
    'heatFluxAvailable',false,'status','未计算','basePressure',NaN, ...
    'basePressureRatio',NaN,'basePressureCoefficient',NaN, ...
    'baseStantonNumber',NaN,'heatFluxRepresentative',NaN, ...
    'heatTransferCoefficientRepresentative',NaN,'heatRecoveryReferenceTemperature',NaN, ...
    'baseEdgeTemperature',NaN,'baseVelocity',NaN,'solvedWallTemperature',NaN, ...
    'shoulderReynolds',NaN,'continuumKnudsen',NaN,'correlationExponentJ',NaN, ...
    'inputStateRelativeResidual',NaN,'heatingCorrelationParameter',NaN, ...
    'totalEnthalpyDifference',NaN,'baseGasTemperatureUsedOnlyForKnScreen',NaN, ...
    'wallTemperatureInput',Twall,'geometry',geometry,'shoulder',shoulder, ...
    'freestream',fs,'pressureMethod','','heatingMethod','','continuumMethod','');
r.reference=struct('report','Lamb and Oberkampf, SAND93-0280, November 1993', ...
    'doi','10.2172/10111104','url','https://www.osti.gov/servlets/purl/10111104', ...
    'pages',[11 21 22 24 26 27 29 30 44 66 67]);
r.heatingScope=['零攻角、无喷流/尾撑、gamma=1.4空气、细长小钝度圆锥；' ...
    '层流9~10度且Re_e,s>=1e5，湍流6~12度；' ...
    '底St实证主要9度钝锥，不能称为独立于建模数据的验证；' ...
    '没有径向分布或面积平均证据，代表底热不能计入全面积解'];
r.geometryVerified=false;
r.eligibleForSurfaceCoverage=false;
r.pressureRequiresMeshAssociation=true;
r.heatFluxRepresentativeOnly=true;
r.radialHeatingDistributionSolved=false;
r.wallTemperatureSolved=false;
r.fieldMetadata=struct;
r.fieldMetadata.wallTemperatureInput=fieldRecord('input','给定壁温边界条件', ...
    isfinite(Twall),'壁温输入，不是本模块求解值');
r.fieldMetadata.boundaryLayerRegime=fieldRecord('input', ...
    '调用者声明的肩部边界层类别；应由其边界层/转捩模块提供',true, ...
    '本底部关联不求解前体转捩，不把类别输入标成求解值');
for key={'basePressure','baseStantonNumber','heatFluxRepresentative', ...
        'heatRecoveryReferenceTemperature','baseEdgeTemperature','baseVelocity', ...
        'solvedWallTemperature'}
    r.fieldMetadata.(key{1})=fieldRecord('unresolved','未求解',false, ...
        '本次尚无适用计算结果');
end
end

function r=fieldRecord(kind,source,valid,scope)
% 模块用途：为每个输出单独保留计算/输入/未解身份与来源。
r=struct('valueKind',kind,'source',source,'valid',logical(valid),'scope',scope);
end
