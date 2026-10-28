function edge = recoverEntropyLayerEdgeState(freestream,station,shockTable,cfg)
%RECOVERENTROPYLAYEREDGESTATE 用激波流管质量坐标匹配边界层外缘熵。
% 模块作用：独立求解 NASA CR-112180 (1972) 式68–69：
%   Psi_BL = rho_e U_e h delta_m，Psi_sh = Psi_BL，delta_m = delta-delta*。
% 对试探 Psi，仅在输入表内线性插值激波穿越角 Gamma；调用同一热物性
% 的前向激波/熵恢复核，得到输入物理静压下的 rho_e、U_e，再闭合质量。
% Gamma 是激波与上游速度夹角，不是壁面转角。总焓由同模型上游推导。
% 本模块不生成激波、物理压力、h、质量厚度或分离判据，不接入默认外流。
%
% station 必需字段：physicalPressure [Pa]、streamlineMetric(h) [m]、
% massThickness(delta_m) [m]、separated (0/1)，以及对应非空来源文字：
% physicalPressureSource、streamlineMetricSource、massThicknessSource、
% separationSource。数值字段须同尺寸或为标量；来源文字可标量广播。
% h 和 delta_m 不得负；零值作为明确的零质量极限单独处理，不用 epsilon。
% shockTable 必需字段：massCoordinate [kg/s] 严格增加；
% shockCrossingAngleDeg 非增加；massCoordinateSource、shockShapeSource、
% upstreamState、upstreamStateSource、thermalPropertyMode。
% upstreamState 必须与 freestream 相同，热物性模式必须与 cfg 相同。
%
% cfg.entropyLayer 可设置 scanSubdivisions(默认16)、relativeTolerance
% (默认1e-9)、maxIterations(默认80)。每个表段内部采样，只有相邻有效
% 样本才形成括区；不外推、不裁剪、不跨越 NaN。数值唯一性表示在报告
% 的扫描分辨率下找到一个候选根，不是数学唯一性证明。多根拒绝输出。
% 结果保留输入/计算来源、残差、表段和实际二分括区；失败场量保持 NaN。

%% 1. 数组、来源和配置合同；结构错误不能悄悄转成物理解
if nargin<4, cfg=struct; end
if ~isstruct(station) || ~isscalar(station) || ~isstruct(shockTable) || ...
        ~isscalar(shockTable) || ~isstruct(cfg) || ~isscalar(cfg)
    error('RocketAeroHeat:EntropyLayer:Input','station、shockTable 和 cfg 必须是标量结构体。');
end
inputNames={'physicalPressure','streamlineMetric','massThickness','separated'};
shape=[1 1];
for j=1:numel(inputNames)
    key=inputNames{j};
    if ~isfield(station,key)
        error('RocketAeroHeat:EntropyLayer:Input','缺少 station.%s。',key);
    end
    if (~isnumeric(station.(key)) && ~islogical(station.(key))) || ...
            ~isreal(station.(key)) || isempty(station.(key))
        error('RocketAeroHeat:EntropyLayer:Input','station.%s 必须为非空实数数组。',key);
    end
    if ~isscalar(station.(key))
        if ~isequal(shape,[1 1]) && ~isequal(shape,size(station.(key)))
            error('RocketAeroHeat:EntropyLayer:Size','站点数组须同尺寸或为标量。');
        end
        shape=size(station.(key));
    end
end
p=expand(station.physicalPressure,shape);
h=expand(station.streamlineMetric,shape);
delta=expand(station.massThickness,shape);
separated=expand(station.separated,shape);
sourceNames={'physicalPressureSource','streamlineMetricSource','massThicknessSource','separationSource'};
sources=struct; sourceValid=true(shape);
for j=1:numel(sourceNames)
    [sources.(sourceNames{j}),goodSource]=stationSource(station,sourceNames{j},shape);
    sourceValid=sourceValid & goodSource;
end
options=struct('scanSubdivisions',16,'relativeTolerance',1e-9,'maxIterations',80);
if isfield(cfg,'entropyLayer')
    if ~isstruct(cfg.entropyLayer) || ~isscalar(cfg.entropyLayer)
        error('RocketAeroHeat:EntropyLayer:Config','cfg.entropyLayer 必须是标量结构体。');
    end
    keys=fieldnames(options);
    for j=1:numel(keys)
        if isfield(cfg.entropyLayer,keys{j}), options.(keys{j})=cfg.entropyLayer.(keys{j}); end
    end
end
validateattributes(options.scanSubdivisions,{'numeric'},{'real','finite','scalar','integer','>=',4,'<=',256});
validateattributes(options.relativeTolerance,{'numeric'},{'real','finite','scalar','>',0,'<',1e-3});
validateattributes(options.maxIterations,{'numeric'},{'real','finite','scalar','integer','positive'});
[mode,thermalScope]=selectAirThermalModel(freestream,cfg);
thermalScope=standaloneAirThermalScope(mode);
edge=initialResult(shape,p,h,delta,separated,sources,mode,thermalScope,options);
inputValid=isfinite(p) & p>0 & isfinite(h) & h>=0 & isfinite(delta) & delta>=0 & ...
    isfinite(separated) & (separated==0 | separated==1);
edge.statusCode(inputValid & ~sourceValid)=2;
usable=inputValid & sourceValid;
edge.statusCode(usable & separated==1)=5;
usable=usable & separated==0;
for j=1:numel(inputNames)
    key=inputNames{j};
    edge.inputMetadata.(key)=struct('valueKind','input','source',{sources.(sourceNames{j})}, ...
        'valid',inputValid & sourceValid,'unit',edge.inputUnits.(key));
end

%% 2. 激波表必须属于同工况、同模型；输入表不允许排序修复或外推
[tableValid,tableMessage,psi,gamma]=checkTable(shockTable,freestream,mode);
edge.tableValidationMessage=tableMessage;
if ~tableValid
    edge.statusCode(usable)=3; edge=finishResult(edge); return;
end
edge.massCoordinateSource=char(shockTable.massCoordinateSource);
edge.shockShapeSource=char(shockTable.shockShapeSource);
edge.upstreamStateSource=char(shockTable.upstreamStateSource);
edge.totalEnthalpySource=['由同热物性上游静焓和完整速度计算；' edge.upstreamStateSource];
if isfield(shockTable,'massCoordinateUnit')
    edge.massCoordinateConvention=char(shockTable.massCoordinateUnit);
end
if isfield(shockTable,'available'), edge.inputTableAvailable=logical(shockTable.available); end
edge.inputTableMassCoordinate=psi; edge.inputTableShockCrossingAngleDeg=gamma;
edge.tableMassCoordinateRange=[psi(1) psi(end)];
if ~any(usable(:)), edge=finishResult(edge); return; end

%% 3. 共用内部扫描：同一激波角的前向跃迁只计算一次
% 先用原始上游静压运行前向核；这里使用 shock 中的总压，并非把该
% 静压误当站点物理压力。各站点之后以各自输入静压重新恢复 rho/U。
scanPsi=psi(1);
for j=1:numel(psi)-1
    part=linspace(psi(j),psi(j+1),options.scanSubdivisions+1)';
    scanPsi=[scanPsi;part(2:end)]; %#ok<AGROW>
end
scanGamma=interp1(psi,gamma,scanPsi,'linear');
kernelCfg=cfg;
if ~isfield(kernelCfg,'externalFlow'), kernelCfg.externalFlow=struct; end
kernelCfg.externalFlow.shockCrossingAngleSource=['输入流管质量/激波角表的表内插值；' edge.shockShapeSource];
kernelCfg.externalFlow.physicalPressureSource='前向扫描用上游静压；站点恢复另用各自输入物理静压';
forward=recoverShockEntropyEdgeState(freestream,scanGamma,freestream.pressure,kernelCfg);
edge.scanSampleCount=numel(scanPsi);
edge.forwardScanUniqueAngleCount=numel(unique(scanGamma));
edge.scanMassCoordinate=scanPsi;
if ~isfield(forward,'upstreamTotalEnthalpy') || ~any(forward.shock.valid(:))
    edge.statusCode(usable)=4; edge=finishResult(edge); return;
end
edge.upstreamTotalEnthalpy=forward.upstreamTotalEnthalpy;
edge.upstreamTotalTemperature=forward.upstreamTotalTemperature;
edge.upstreamTotalPressure=forward.upstreamTotalPressure;

%% 4. 求解标量质量方程；零质量、表端点、多根和无效域分别处理
for index=reshape(find(usable),1,[])
    C=h(index)*delta(index);
    scanState=recoverCachedForward(freestream,forward,p(index),mode);
    valid=scanState.valid & forward.shock.valid;
    boundaryMass=C*scanState.density.*scanState.speed;
    residual=scanPsi-boundaryMass;
    residual(~valid)=NaN;
    edge.validScanSampleCount(index)=nnz(valid);
    edge.invalidScanIntervalCount(index)=nnz(~(valid(1:end-1) & valid(2:end)));
    if ~any(valid)
        edge.statusCode(index)=4; continue;
    end
    % h=0 或质量厚度=0 是输入给出的退化极限；Psi=0 必须真在表内。
    if C==0
        if psi(1)~=0 || ~valid(1)
            edge.statusCode(index)=6; continue;
        end
        candidates=[0 0];
        edge.zeroMassLimit(index)=true;
    else
        scale=max(abs(scanPsi),abs(boundaryMass));
        exact=valid & abs(residual)<=64*eps.*max(scale,realmin);
        candidates=[scanPsi(exact) scanPsi(exact)];
        crossings=find(valid(1:end-1) & valid(2:end) & ~exact(1:end-1) & ~exact(2:end) & ...
            sign(residual(1:end-1))~=sign(residual(2:end)));
        candidates=[candidates;scanPsi(crossings) scanPsi(crossings+1)]; %#ok<AGROW>
        if ~isempty(candidates), candidates=sortrows(candidates,1); end
    end
    edge.rootCountDetected(index)=size(candidates,1);
    edge.candidateMassBrackets{index}=candidates;
    if isempty(candidates)
        edge.statusCode(index)=6; continue;
    elseif size(candidates,1)>1
        edge.statusCode(index)=7; continue;
    end
    edge.numericalUniquenessChecked(index)=true;
    edge.scanBracketLower(index)=candidates(1,1);
    edge.scanBracketUpper(index)=candidates(1,2);
    lo=candidates(1,1); hi=candidates(1,2);
    [left,fl]=evaluateMass(lo,p(index),C,freestream,psi,gamma,kernelCfg);
    if ~left.algebraicallyValid, edge.statusCode(index)=4; continue; end
    root=lo; result=left; converged=massConverged(root,fl,options.relativeTolerance);
    if ~converged && hi>lo
        [right,fr]=evaluateMass(hi,p(index),C,freestream,psi,gamma,kernelCfg);
        if ~right.algebraicallyValid || sign(fl)==sign(fr)
            edge.statusCode(index)=8; continue;
        end
        for iteration=1:options.maxIterations
            root=lo+(hi-lo)/2;
            [result,fm]=evaluateMass(root,p(index),C,freestream,psi,gamma,kernelCfg);
            edge.iterations(index)=iteration;
            if ~result.algebraicallyValid
                % 不把 NaN 当成带符号数，也不跨过物性空洞继续二分。
                edge.statusCode(index)=4; break;
            end
            if massConverged(root,fm,options.relativeTolerance)
                converged=true; break;
            end
            if sign(fm)==sign(fl), lo=root; fl=fm; else, hi=root; end
            if root==lo && lo==hi, break; end
        end
    end
    edge.finalBracketLower(index)=lo; edge.finalBracketUpper(index)=hi;
    if ~converged
        if edge.statusCode(index)~=4, edge.statusCode(index)=8; end
        continue;
    end
    matchedMass=C*result.density*result.speed;
    difference=root-matchedMass;
    if ~massConverged(root,difference,options.relativeTolerance)
        edge.statusCode(index)=8; continue;
    end
    outputNames={'pressure','temperature','density','speed','Mach','speedOfSound', ...
        'totalPressure','totalTemperature','entropyIncrease','energyResidual','entropyResidual'};
    for j=1:numel(outputNames), edge.(outputNames{j})(index)=result.(outputNames{j}); end
    edge.massCoordinate(index)=root;
    edge.boundaryLayerMassCoordinate(index)=matchedMass;
    edge.shockCrossingAngleDeg(index)=interp1(psi,gamma,root,'linear');
    edge.massResidual(index)=difference;
    edge.massResidualNormalized(index)=difference/max([abs(root),abs(matchedMass),realmin]);
    segment=find(psi<=root,1,'last');
    if segment==numel(psi), segment=segment-1; end
    edge.tableSegmentIndex(index)=segment;
    edge.tableBracketLower(index)=psi(segment);
    edge.tableBracketUpper(index)=psi(segment+1);
    edge.statusCode(index)=0;
    edge.algebraicallyValid(index)=true;
    edge.massMatchingSolved(index)=true;
    edge.stagnationLimit(index)=result.speed==0;
    edge.normalKernelMethod{index}=result.shock.normalKernelMethod{1};
end
edge=finishResult(edge);
end

function [value,residual]=evaluateMass(psiValue,p,C,fs,psi,gamma,cfg)
% 真正迭代新穿越角；物理压力的来源与扫描辅助压力分别保留。
angle=interp1(psi,gamma,psiValue,'linear');
cfg.externalFlow.physicalPressureSource='质量匹配站点输入物理静压；本内核未求压力分布';
value=recoverShockEntropyEdgeState(fs,angle,p,cfg);
residual=psiValue-C*value.density*value.speed;
end

function yes=massConverged(psiValue,residual,tolerance)
% 零端点不挪动；一般根采用相对质量残差，不使用任意绝对 epsilon。
scale=max([abs(psiValue),abs(psiValue-residual),realmin]);
yes=isfinite(residual) && abs(residual)<=tolerance*scale;
end

function state=recoverCachedForward(fs,forward,p,mode)
% 复用所有共用扫描角的激波后总压；仅重新反解该站点物理静压。
p0=forward.shock.totalPressure(:); T0=forward.upstreamTotalTemperature;
if strcmp(mode,'thermallyPerfect')
    state=recoverThermallyPerfectEdgeState(repmat(p,size(p0)),p0,T0,fs.R);
else
    valid=isfinite(p0) & p>0 & p<=p0;
    cp=fs.gamma*fs.R/(fs.gamma-1);
    T=nan(size(p0)); U=T;
    T(valid)=T0*(p./p0(valid)).^((fs.gamma-1)/fs.gamma);
    U(valid)=sqrt(2*cp*(T0-T(valid)));
    state=struct('valid',valid,'temperature',T,'density',p./(fs.R*T),'speed',U);
    state.density(~valid)=NaN; state.speed(~valid)=NaN;
end
% NASA7 总压等于物理静压时使用准确零速度，避免逆解舍入噪声。
atStagnation=state.valid & p==p0;
state.temperature(atStagnation)=T0;
state.density(atStagnation)=p/(fs.R*T0);
state.speed(atStagnation)=0;
end

function [valid,message,psi,gamma]=checkTable(table,fs,mode)
% 同一工况签名包含速度及密度，不能仅以名义 Mach 相同判定一致。
valid=false; psi=[]; gamma=[]; message='激波表字段或来源缺失';
% 生成器明确拒绝的表不能绕过 available 状态被当作已生成结果消费。
if isfield(table,'available')
    flag=table.available;
    if (~isnumeric(flag) && ~islogical(flag)) || ~isreal(flag) || ...
            ~isscalar(flag) || ~isfinite(flag) || flag~=1
        message='输入激波表显式标为不可用'; return;
    end
end
if isfield(table,'massCoordinateUnit') && ~sourceText(table.massCoordinateUnit)
    message='输入激波质量坐标单位/流管归一化约定无效'; return;
end
keys={'massCoordinate','shockCrossingAngleDeg','massCoordinateSource','shockShapeSource', ...
    'upstreamState','upstreamStateSource','thermalPropertyMode'};
if ~all(isfield(table,keys)), return; end
for key={'massCoordinateSource','shockShapeSource','upstreamStateSource','thermalPropertyMode'}
    if ~sourceText(table.(key{1})), message='激波表来源/热物性模式必须为非空标量文字'; return; end
end
if ~strcmp(char(table.thermalPropertyMode),mode)
    message='激波表与恢复内核的热物性模式不同'; return;
end
if ~isnumeric(table.massCoordinate) || ~isreal(table.massCoordinate) || ...
        ~isvector(table.massCoordinate) || ~isnumeric(table.shockCrossingAngleDeg) || ...
        ~isreal(table.shockCrossingAngleDeg) || ~isvector(table.shockCrossingAngleDeg)
    message='质量坐标和激波角须为实数向量'; return;
end
psi=double(table.massCoordinate(:)); gamma=double(table.shockCrossingAngleDeg(:));
if numel(psi)<2 || numel(psi)~=numel(gamma) || ~all(isfinite(psi)) || ...
        ~all(isfinite(gamma)) || psi(1)<0 || any(diff(psi)<=0) || ...
        any(diff(gamma)>0) || any(gamma<=0 | gamma>90)
    message='质量坐标须非负严格增加，激波角须位于(0,90]且非增加；不排序、不裁剪'; return;
end
if ~isstruct(table.upstreamState) || ~isscalar(table.upstreamState)
    message='激波表缺少同工况上游结构'; return;
end
stateKeys={'Mach','gamma','R','pressure','temperature','density','speed'};
% 几何表在零攻角生成时，旋转工况不能沿用相同标量速度的旧激波表。
for key={'alpha','beta'}
    if isfield(fs,key{1}) || isfield(table.upstreamState,key{1})
        stateKeys{end+1}=key{1}; %#ok<AGROW>
    end
end
for key=stateKeys
    name=key{1};
    if ~isfield(fs,name) || ~isfield(table.upstreamState,name)
        message='激波表或调用上游缺少热力状态字段'; return;
    end
    a=fs.(name); b=table.upstreamState.(name);
    isAngle=ismember(name,{'alpha','beta'});
    if ~isnumeric(a) || ~isreal(a) || ~isscalar(a) || ~isfinite(a) || (~isAngle && a<=0) || ...
            ~isnumeric(b) || ~isreal(b) || ~isscalar(b) || ~isfinite(b) || (~isAngle && b<=0) || ...
            abs(a-b)>1e-10*max(abs(a),abs(b))
        message=['激波表上游与调用上游不一致：' name]; return;
    end
end
valid=true; message='输入表结构、单调性、来源及同工况/热物性签名已通过';
end

function result=initialResult(shape,p,h,delta,separated,sources,mode,thermalScope,options)
% 每个计算量先置 NaN；原始输入单独保留，不能用输入冒充成功计算。
blank=nan(shape);
result=struct('available',false,'algebraicallyValid',false(shape), ...
    'massMatchingSolved',false(shape),'geometryVerified',false, ...
    'defaultFlowIntegrated',false,'statusCode',ones(shape,'uint8'), ...
    'statusNames',{{'已计算（独立表内质量匹配）','站点数值输入无效', ...
    '站点来源缺失或无效','激波表/来源/同工况/热物性合同无效', ...
    '前向激波、物理静压或热物性域无效','分离区不适用', ...
    '没有检测到表内有效质量匹配根','检测到多个根，不能选择唯一外缘', ...
    '括区或迭代未收敛'}},'thermalPropertyMode',mode,'thermalScope',thermalScope, ...
    'physicalPressureInput',p,'streamlineMetricInput',h,'massThicknessInput',delta, ...
    'separatedInput',separated,'inputSources',sources,'inputMetadata',struct, ...
    'inputUnits',struct('physicalPressure','Pa','streamlineMetric','m', ...
    'massThickness','m','separated','logical'), ...
    'massCoordinateUnit','kg/s','scanSubdivisionsPerTableSegment',options.scanSubdivisions, ...
    'massCoordinateConvention','kg/s；输入质量坐标和h必须采用同一流管归一化', ...
    'inputTableAvailable',NaN, ...
    'relativeMassTolerance',options.relativeTolerance,'maximumIterations',options.maxIterations, ...
    'scanSampleCount',0,'forwardScanUniqueAngleCount',0, ...
    'scanMassCoordinate',[],'tableMassCoordinateRange',[NaN NaN], ...
    'inputTableMassCoordinate',[],'inputTableShockCrossingAngleDeg',[], ...
    'rootCountDetected',zeros(shape),'numericalUniquenessChecked',false(shape), ...
    'uniquenessMethod','逐表段内部均匀采样 + 有效相邻样本符号括区；数值检查，非数学唯一性证明', ...
    'validScanSampleCount',zeros(shape),'invalidScanIntervalCount',zeros(shape), ...
    'iterations',zeros(shape),'candidateMassBrackets',{cell(shape)}, ...
    'normalKernelMethod',{cell(shape)},'zeroMassLimit',false(shape),'stagnationLimit',false(shape), ...
    'massCoordinateSource','','shockShapeSource','','upstreamStateSource','', ...
    'totalEnthalpySource','','upstreamTotalEnthalpy',NaN,'upstreamTotalTemperature',NaN, ...
    'upstreamTotalPressure',NaN,'tableValidationMessage','尚未检查输入表', ...
    'method',['NASA CR-112180 式68–69标量质量匹配 + 表内激波角插值 + ' ...
    '前向法向跃迁/切向连续/同总焓熵恢复；' mode], ...
    'scope',['仅独立附着流站点；物理压力、流管度量、质量厚度和激波表均为输入；' ...
    '不验证自算几何、表的物理准确性、分离再附、热流或真实鼻部/后体匹配；' ...
    '不外推且数值唯一性受报告的扫描分辨率限制'], ...
    'sourceURL','https://ntrs.nasa.gov/api/citations/19730005231/downloads/19730005231.pdf', ...
    'status','尚未完成独立质量匹配');
fields={'pressure','temperature','density','speed','Mach','speedOfSound','totalPressure', ...
    'totalTemperature','entropyIncrease','energyResidual','entropyResidual','massCoordinate', ...
    'boundaryLayerMassCoordinate','shockCrossingAngleDeg','massResidual','massResidualNormalized', ...
    'tableSegmentIndex','tableBracketLower','tableBracketUpper','scanBracketLower','scanBracketUpper', ...
    'finalBracketLower','finalBracketUpper'};
result.fieldMetadata=struct;
for j=1:numel(fields)
    key=fields{j}; result.(key)=blank;
    metadata=struct('valueKind','computed','source',result.method,'scope',result.scope, ...
        'valid',false(shape),'statusCode',result.statusCode);
    if strcmp(key,'pressure')
        metadata.valueKind='input'; metadata.source=sources.physicalPressureSource;
    end
    result.fieldMetadata.(key)=metadata;
end
end

function result=finishResult(result)
% 所有早退和常规出口统一同步逐变量状态，不遗留初始化有效标记。
fields=fieldnames(result.fieldMetadata);
for j=1:numel(fields)
    result.fieldMetadata.(fields{j}).valid=result.algebraicallyValid;
    result.fieldMetadata.(fields{j}).statusCode=result.statusCode;
end
result.entropySource=['求得穿越质量坐标对应的表内 Gamma；' result.shockShapeSource ...
    '；前向激波熵及同总焓恢复，不沿后体人为平滑'];
result.available=any(result.algebraicallyValid(:));
result.status=sprintf('独立质量匹配 %d/%d 个站点；数值扫描不是数学唯一性证明；未接入默认外流', ...
    nnz(result.algebraicallyValid),numel(result.algebraicallyValid));
end

function value=expand(input,shape)
value=double(input); if isscalar(value), value=repmat(value,shape); end
end

function yes=sourceText(value)
yes=(ischar(value) && (isrow(value) || isempty(value))) || (isstring(value) && isscalar(value));
if yes, yes=~ismissing(string(value)) && strlength(strtrim(string(value)))>0; end
end

function [values,valid]=stationSource(station,key,shape)
values=repmat({''},shape); valid=false(shape);
if ~isfield(station,key), return; end
raw=station.(key);
if sourceText(raw), values(:)={char(raw)}; valid(:)=true; return; end
if isstring(raw), raw=cellstr(raw); end
if ~iscell(raw) || ~isequal(size(raw),shape), return; end
for j=1:numel(raw)
    if sourceText(raw{j}), values{j}=char(raw{j}); valid(j)=true; end
end
end
