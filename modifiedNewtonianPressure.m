function mesh = modifiedNewtonianPressure(mesh,freestream,cfg)
%MODIFIEDNEWTONIANPRESSURE 计算修正牛顿面元压力，结果存入 mesh.pressure。
% 输入：最终网格及 localGeometry、当前来流、可选 cfg。
% 迎风公式：Cp=CpMax*max(mu,0)^2，p=pinf+qinf*Cp。
% 背风及切向默认 Cp=0，即自由来流静压；这不是分离流/尾迹压力模型。
% 可设置 cfg.pressure.leewardCp 和 cfg.pressure.aftBaseCp 为常数。
% aftBaseCp 仅覆盖非迎风尾部平封盖；不覆盖反向流中的迎风尾盖。
% cfg.basePressure.model='sand93Laminar' 可选择经真实几何/肩部薄层审核
% 的尖锥底压关联；不支持的尾底保持NaN，不暗中退回给定Cp。
% 方法主要用于高超声速迎风面工程估算；M>=1 只是公式定义域。
% 不处理凹形表面遮挡、激波干扰、粘性、真实气体或气动热。

if nargin < 3, cfg = struct; end
%% 1. 计算参考上限，重新分类以避免使用上一次攻角的分类结果
[CpMax,reference] = calcCpMax(freestream);
mesh = classifyAerothermalRegion(mesh,freestream,cfg);
r = mesh.aerothermalRegion;
leewardCp = 0; aftBaseCp = 0;
if isfield(cfg,'pressure')
    if isfield(cfg.pressure,'leewardCp'), leewardCp = cfg.pressure.leewardCp; end
    if isfield(cfg.pressure,'aftBaseCp'), aftBaseCp = cfg.pressure.aftBaseCp; end
end
validateattributes(leewardCp,{'numeric'},{'real','finite','scalar'});
validateattributes(aftBaseCp,{'numeric'},{'real','finite','scalar'});

%% 2. 使用原始 mu 而不是分类阈值，保留略微迎风面的非零压力增量
mu = r.mu; windward = mu > 0;
Cp = leewardCp*ones(size(mu));
Cp(windward) = CpMax*mu(windward).^2;
base = r.isAftCap & ~windward;
Cp(base) = aftBaseCp;
qinf = reference.dynamicPressure;
pinf = freestream.pressure;
pressure = pinf+qinf*Cp;
if any(~isfinite(pressure) | pressure < 0)
    error('RocketAeroHeat:Pressure:Nonphysical','指定压力系数产生负绝对压力或非有限值。');
end

%% 2a. 独立底区模型：先审核再写入当前压力，失败区域保留未解
% 默认仍采用明确给定的底Cp。选择关联时，不把给定边界当作关联解；
% 所有实际底面值、来源和拒绝原因跟随本次工况保存，供积分/网页使用。
baseModel='prescribed'; baseRegion=[];
if isfield(cfg,'basePressure') && isfield(cfg.basePressure,'model')
    baseModel=char(cfg.basePressure.model);
end
if ~ismember(baseModel,{'prescribed','sand93Laminar'})
    error('RocketAeroHeat:Pressure:BaseModel','未知底压模型：%s。',baseModel);
end
fieldStatus=ones(size(mu),'uint8'); fieldStatus(~windward)=uint8(2);
faceMethod=repmat({'修正牛顿法（定比热完全气体）'},numel(mu),1);
faceScope=repmat({'迎风高超声速工程近似；不含遮挡、激波干扰或粘性压力'},numel(mu),1);
faceMethod(~windward & ~base)={'给定背风Cp边界条件'};
faceScope(~windward & ~base)={'给定背风边界；不是分离流压力求解'};
faceMethod(base)={'给定尾底Cp边界条件'};
faceScope(base)={'给定非迎风尾盖边界；不是尾迹求解'};
computedBase=false(size(base)); prescribedBase=base;
if strcmp(baseModel,'sand93Laminar')
    baseRegion=solveMeshBasePressure(mesh,freestream,cfg);
    if ~baseRegion.geometry.flatBaseVerified
        error('RocketAeroHeat:Pressure:BaseTarget', ...
            ['当前尾底未通过完整平圆底目标识别：%s。不能安全定位全部拒绝面元，' ...
             '所选底压模型未应用；请保留明确给定边界或更换有支持的几何。'], ...
            baseRegion.geometry.status);
    end
    target=base | (baseRegion.baseFace & ~windward);
    Cp(target)=NaN;pressure(target)=NaN;fieldStatus(target)=uint8(0);
    prescribedBase(target)=false;
    faceMethod(target)={['底压未解：' baseRegion.status]};
    faceScope(target)={'所选SAND93-0280关联未通过真实几何、自由流、操作声明或肩部薄层审核；没有给定Cp回退'};
    computedBase=baseRegion.computedFace & target;
    Cp(computedBase)=baseRegion.Cp(computedBase);
    pressure(computedBase)=baseRegion.pressure(computedBase);
    fieldStatus(computedBase)=uint8(1);
    faceMethod(computedBase)=baseRegion.faceMethod(computedBase);
    faceScope(computedBase)=baseRegion.faceScope(computedBase);
end

%% 3. 保存压力、近似来源与计算时的来流，供 GUI 和后续积分读取
result = struct;
result.model = '修正牛顿法（定比热完全气体）';
result.Cp = Cp; result.CpMax = CpMax;
result.staticPressure = pressure; result.pressureIncrement = qinf*Cp;
result.freestreamPressure = pinf; result.dynamicPressure = qinf;
result.reference = reference; result.freestream = freestream;
result.mu = mu; result.isNewtonianFace = windward;
% 保存压力所对应的几何，供整体载荷模块识别过期压力。
result.geometry = struct('normal',mesh.localGeometry.normal, ...
    'center',mesh.localGeometry.center,'area',mesh.localGeometry.area);
result.isPrescribedBase = prescribedBase; result.isPrescribedLeeward = ~windward & ~base & fieldStatus==2;
result.isComputedBase=computedBase; result.baseModel=baseModel;
result.baseRegion=baseRegion; result.fieldStatus=fieldStatus;
result.faceMethod=faceMethod; result.faceScope=faceScope;
result.missingPressureFace=fieldStatus==0;
result.leewardCp = leewardCp; result.aftBaseCp = aftBaseCp;
result.summary = struct('minimumCp',min(Cp),'maximumCp',max(Cp), ...
    'minimumPressure',min(pressure),'maximumPressure',max(pressure));
mesh.pressure = result;
% 旧载荷不能继续搭配新的压力；调用方需重新积分。
if isfield(mesh,'aerodynamicLoads'), mesh=rmfield(mesh,'aerodynamicLoads'); end
if isfield(mesh,'stagnation'), mesh=rmfield(mesh,'stagnation'); end
if isfield(mesh,'boundaryLayerEdge'), mesh=rmfield(mesh,'boundaryLayerEdge'); end
if isfield(mesh,'boundaryLayer'), mesh=rmfield(mesh,'boundaryLayer'); end
if isfield(mesh,'heating'), mesh=rmfield(mesh,'heating'); end
if isfield(mesh,'heatingCoverage'), mesh=rmfield(mesh,'heatingCoverage'); end
if isfield(mesh,'stagnationHeating'), mesh=rmfield(mesh,'stagnationHeating'); end
if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
if isfield(mesh,'surfaceStreamlines'), mesh=rmfield(mesh,'surfaceStreamlines'); end
if isfield(mesh,'axisymmetricMarch'), mesh=rmfield(mesh,'axisymmetricMarch'); end
if isfield(mesh,'externalFlow')
    saved=mesh.externalFlow.freestream;
    same=true;
    for key={'Mach','alpha','beta','altitude','gamma','R', ...
            'pressure','temperature','density','speed','totalPressure','totalTemperature'}
        name=key{1};
        same=same && isfield(saved,name) && ...
            isequal(saved.(name),freestream.(name));
    end
    if ~same, mesh=rmfield(mesh,'externalFlow'); end
end
end
