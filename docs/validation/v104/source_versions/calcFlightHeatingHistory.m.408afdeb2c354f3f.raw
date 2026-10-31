function result = calcFlightHeatingHistory(mesh,flightHistory,cfg)
%CALCFLIGHTHEATINGHISTORY 在同一网格上逐站计算内置外流和气动热边界条件。
% 模块作用：用户只给时间、海拔、Mach、攻角和侧滑角；每个时刻都从
% 标准大气、内置激波/膨胀、三维准流线、边界层、转捩算到面元 h 与 Taw。
% 输出供 calcWallThermalHistory 使用；绝不把前一站的外流误当当前状态。
% 用户时间站之间按明确的最大步长线性插值飞行条件；每个新站重算
% 定常外流，壁温在相邻计算站之间仍按右端热边界条件推进。
% 这不等于非定常激波或连续轨迹求解。

%% 1. 校验飞行历程和最终网格
if isfield(cfg,'heating') && isfield(cfg.heating,'noseModel') && ...
        strcmp(cfg.heating.noseModel,'sphericalLocalSimilarity')
    error('RocketAeroHeat:FlightHistory:SingleConditionModel', ...
        '球鼻局部相似模型仅验证单工况；暂不能用于飞行历程热边界。');
end
fields={'time','altitude','Mach','alpha','beta'};
for k=1:numel(fields)
    name=fields{k};
    if ~isfield(flightHistory,name)
        error('RocketAeroHeat:FlightHistory:Field','缺少 %s。',name);
    end
    validateattributes(flightHistory.(name),{'numeric'}, ...
        {'real','finite','vector','nonempty'},mfilename,name);
end
t=double(flightHistory.time(:)'); nt=numel(t);
if nt<2 || t(1)~=0 || any(diff(t)<=0)
    error('RocketAeroHeat:FlightHistory:Time', ...
        '时间至少有两个站点，必须从 0 开始并严格递增。');
end
for k=2:numel(fields)
    if numel(flightHistory.(fields{k}))~=nt
        error('RocketAeroHeat:FlightHistory:Length', ...
            '%s 的长度须与时间相同。',fields{k});
    end
end
mach=double(flightHistory.Mach(:)');
if any(mach<=1)
    error('RocketAeroHeat:FlightHistory:Mach', ...
        '当前内置激波/膨胀链路要求各站 Mach>1。');
end
if ~isfield(cfg,'flight') || ~isfield(cfg.flight,'Twall')
    error('RocketAeroHeat:FlightHistory:Config','缺少初始壁温 cfg.flight.Twall。');
end
if ~isfield(mesh,'localGeometry') || ~isfield(mesh,'surfaceCoord')
    error('RocketAeroHeat:FlightHistory:Mesh', ...
        '需要完成几何分析和表面坐标的最终网格。');
end
inputFlight=flightHistory;
maxTimeStep=1;
if isfield(cfg,'flightHistory') && ...
        isfield(cfg.flightHistory,'maxTimeStep')
    maxTimeStep=cfg.flightHistory.maxTimeStep;
end
validateattributes(maxTimeStep,{'numeric'}, ...
    {'real','finite','scalar','positive'});
nf=size(mesh.F,1);
estimatedStations=1+sum(max(1,ceil(diff(t)/maxTimeStep)));
if ~isfinite(estimatedStations) || nf*estimatedStations>1e6
    error('RocketAeroHeat:FlightHistory:Memory', ...
        '面元数乘加密后的计算站数超过 1e6，请增大最大计算间隔。');
end
[flightHistory,stationRefinement]= ...
    refineFlightStations(flightHistory,maxTimeStep);
t=flightHistory.time;
nt=numel(t);
mach=flightHistory.Mach;

%% 2. 逐站独立求解外流和热流边界条件
h=nan(nf,nt); Taw=h; eligible=false(nf,nt);
methodID=zeros(nf,nt,'uint8');
Te=h; p=h; exponent=h; nu0=h; cp0=h;
pointH=nan(1,nt); pointTaw=pointH; pointValid=false(1,nt);
pointPosition=nan(nt,3);
validCounts=zeros(1,nt); noseCounts=validCounts;
detachedCounts=validCounts;
axialCounts=validCounts;
axialApplied=false(1,nt);
axialStatus=cell(1,nt);
areaCoverage=zeros(1,nt);
reasonArea=[];
reasonCode=zeros(nf,nt,'uint8');
knudsenByStation=nan(1,nt);
rarefiedFacesByStation=zeros(1,nt);
unknownThermalFacesByStation=zeros(1,nt);
upstreamUnknownThermalFacesByStation=zeros(1,nt);
freestreams=cell(1,nt);
reusedSteadyState=false(1,nt);
flowStations=cell(1,nt);
initial=cfg.flight.Twall;
for j=1:nt
    condition=cfg.flight;
    condition.altitude=double(flightHistory.altitude(j));
    condition.Mach=mach(j);
    condition.alpha=double(flightHistory.alpha(j));
    condition.beta=double(flightHistory.beta(j));
    condition.Twall=initial;
    stationCondition=[condition.altitude,condition.Mach, ...
        condition.alpha,condition.beta];
    if j==1 || ~isequal(stationCondition,previousStationCondition)
        % 连续两个计算站工况完全相同时，重用同一份定常外流与热流。
        % 壁温积分仍按两个独立时间站执行，不跳过时间步。
        fs=calcFreestream(condition);
        current=calcFlightStationState(mesh,condition,cfg);
        flowSnapshot=captureAerothermalFlowSnapshot(current);
    else
        reusedSteadyState(j)=true;
    end
    previousStationCondition=stationCondition;
    flowStations{j}=flowSnapshot;
    axialApplied(j)=current.axisymmetricMarch.applied;
    axialStatus{j}=current.axisymmetricMarch.status;
    if axialApplied(j)
        axialCounts(j)=current.axisymmetricMarch.summary.marchedFaces;
    end
    heat=current.heating; edge=current.boundaryLayerEdge;
    if j==1
        airMode=heat.thermalPropertyMethod;
        airScope=heat.thermalPropertyScope;
    end
    h(:,j)=heat.heatTransferCoefficient;
    Taw(:,j)=heat.recoveryTemperature;
    eligible(:,j)=heat.eligibleFace;
    methodID(:,j)=heat.faceMethodID;
    validCounts(j)=heat.summary.calculatedFaces;
    noseCounts(j)=heat.summary.noseFaces;
    detachedCounts(j)=heat.summary.detachedClosureHeatFaces;
    coverage=current.heatingCoverage;
    areaCoverage(j)=coverage.summary.computedAreaFraction;
    if j==1, reasonArea=zeros(numel(coverage.reasonNames),nt); end
    reasonArea(:,j)=coverage.areaByReason(:);
    reasonCode(:,j)=coverage.reasonCode;
    knudsenByStation(j)=current.continuumRegime.freestreamKnudsen;
    rarefiedFacesByStation(j)= ...
        current.continuumRegime.summary.directKnudsenUnresolvedFaces;
    % 兼容已有历程时也按真实身份统计；未知占位不写成稀薄流面数。
    unknownThermalFacesByStation(j)= ...
        current.continuumRegime.summary.directStateUnresolvedFaces;
    upstreamUnknownThermalFacesByStation(j)= ...
        current.continuumRegime.summary.downstreamStateUnresolvedFaces;
    freestreams{j}=fs;
    if current.stagnationHeating.available
        pointValid(j)=true;
        pointH(j)=current.stagnationHeating.heatTransferCoefficient;
        pointTaw(j)=current.stagnationHeating.adiabaticWallTemperature;
        pointPosition(j,:)=current.stagnationHeating.position(:)';
    end

    % Eckert 参考性质反馈以初始壁温为基准；鼻端余弦延拓保持固定 h。
    localExponent=nan(nf,1);
    localExponent(heat.laminarFace)=0.5;
    localExponent(heat.turbulentFace)=0.2;
    localExponent(heat.noseFace)=0;
    localTe=edge.temperature(:); localP=edge.pressure(:);
    localTe(heat.noseFace)=Taw(heat.noseFace,j);
    localP(heat.noseFace)=fs.pressure;
    reference=localTe+0.5*(initial-localTe)+ ...
        0.22*(Taw(:,j)-localTe);
    mu=1.458e-6*reference.^(3/2)./(reference+110.4);
    localNu=mu.*fs.R.*reference./localP;
    active=eligible(:,j);
    if any(~isfinite(localNu(active)) | localNu(active)<=0 | ...
            ~isfinite(localExponent(active)))
        error('RocketAeroHeat:FlightHistory:Reference', ...
            '第 %d 站的参考性质或算法类别无效。',j);
    end
    Te(:,j)=localTe; p(:,j)=localP;
    exponent(:,j)=localExponent; nu0(:,j)=localNu;
    if strcmp(airMode,'thermallyPerfect')
        localCp=airEnthalpySecantCp(Taw(:,j),initial);
        if any(~isfinite(localCp(active)) | localCp(active)<=0)
            error('RocketAeroHeat:FlightHistory:AirRange', ...
                '第 %d 站的壁温或恢复温度超出干空气焓模型适用域。',j);
        end
        cp0(:,j)=localCp;
    end
end

%% 3. 保存每站算法来源、适用性和壁面边界条件
result=struct;
if strcmp(airMode,'thermallyPerfect')
    airLabel='干空气 NASA7/低温振动冻结焓恢复';
else
    airLabel='定比热恢复温度';
end
result.method=[stationRefinement.method ' + 逐站内置局部切锥/切楔/膨胀外流 + ' ...
    '零攻角轴向波系推进 + 三维准流线 + ' ...
    '受限前体脱体闭合 + Thwaites/经验转捩 + Eckert/Reynolds 热流 + ' ...
    'Sutton-Graves 球鼻偏置驻点余弦延拓 + Kn 连续介质筛选 + ' airLabel];
result.scope=['同一刚性网格、各计算站定常工程近似；用户站之间飞行条件' ...
    '分段线性插值；壁温积分采用相邻计算站右端热边界条件；' ...
    '鼻端仅在已验证球形几何与合成入射角 <=30° 时低阶估算；' ...
    '前体脱体闭合仅限局部切楔超限、前 10% 长度与转角 <=60°；' ...
    '轴向波系推进仅限零攻角且能识别近轴对称直鼻锥；' ...
    '疑似层流分离及其下游因无再附模型而留空；' ...
    'Kn 超阈值的滑移/过渡流不套用连续介质热流；' ...
    '此处 h/Taw 和覆盖率采用统一初始壁温，瞬态壁温耦合由壁温求解器另行记录；' ...
    '无非定常激波、真实气体或烧蚀；' airScope];
result.time=t;
result.flight=flightHistory;
result.inputFlight=inputFlight;
result.stationRefinement=stationRefinement;
result.reusedSteadyState=reusedSteadyState;
result.flowStations=flowStations;
result.initialWallTemperature=initial;
% 保存原配置和精确网格顺序，用于按需重建耦合时间站；不缓存全部
% 外流大结构，防止多站复制几何/路径/派生量造成不必要的内存增长。
result.couplingInput=struct('configuration',cfg,'vertices',mesh.V,'faces',mesh.F);
result.heatTransferCoefficient=h;
result.recoveryTemperature=Taw;
result.eligibleFace=eligible;
result.faceMethodID=methodID;
result.axisymmetricMarchAppliedByStation=axialApplied;
result.axisymmetricMarchStatusByStation=axialStatus;
result.faceReasonCode=reasonCode;
result.reasonNames=coverage.reasonNames;
result.referenceLaw=struct('edgeTemperature',Te,'edgePressure',p, ...
    'exponent',exponent,'referenceViscosity0',nu0, ...
    'gasConstant',freestreams{1}.R);
if strcmp(airMode,'thermallyPerfect')
    result.referenceLaw.thermalPropertyMode=airMode;
    result.referenceLaw.referenceCp0=cp0;
end
result.stagnation=struct('available',pointValid, ...
    'heatTransferCoefficient',pointH,'recoveryTemperature',pointTaw, ...
    'position',pointPosition);
result.freestream=freestreams;
result.freestreamKnudsenByStation=knudsenByStation;
result.rarefiedFacesByStation=rarefiedFacesByStation;
result.unknownThermalFacesByStation=unknownThermalFacesByStation;
result.upstreamUnknownThermalFacesByStation=upstreamUnknownThermalFacesByStation;
result.continuumMethod=current.continuumRegime.method;
result.continuumScope=current.continuumRegime.scope;
result.summary=struct('calculatedFacesByStation',validCounts, ...
    'inputStations',stationRefinement.inputStationCount, ...
    'evaluationStations',stationRefinement.evaluationStationCount, ...
    'reusedSteadyStates',nnz(reusedSteadyState), ...
    'noseFacesByStation',noseCounts, ...
    'detachedClosureHeatFacesByStation',detachedCounts, ...
    'axialMarchFacesByStation',axialCounts, ...
    'calculatedAreaFractionByStation',areaCoverage, ...
    'reasonAreaByStation',reasonArea, ...
    'validAcrossHistory',nnz(all(eligible,2)), ...
    'validAcrossHistoryAreaFraction', ...
    sum(mesh.localGeometry.area(all(eligible,2)))/ ...
        sum(mesh.localGeometry.area));
end
