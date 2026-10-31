function state = aeroHeatWorkflowState(action,state,payload)
%AEROHEATWORKFLOWSTATE 五步骤工作台的集中状态控制器（不执行物理计算）。
% 模块一：保存 SI 几何、已确认朝向、网格、配置与计算快照；修改上游
% 输入时使下游失效。页面浏览不计算、不更改输入版本。
% 模块二：显式求解时冻结输入版本和求解编号，只接受对应当前版本的
% 完整快照；失败或过期结果不能冒充当前计算完成。
% 模块三：生成侧边栏前置条件及完成状态，供同一宿主界面使用。
%
% state = aeroHeatWorkflowState(action,state,payload)
% init:             payload.cfg（可省略，使用空配置）
% setModel:         payload.rawMesh；可选 orientedMesh、cfg、modelSource
% invalidateModel:  文件/尺度编辑后清模型及下游；可选 payload.cfg
% setOrientation:   payload.orientedMesh；可选 cfg
% invalidateOrientation: 朝向预览实际变更，保留 raw；可选 payload.cfg
% setMesh:          payload.mesh；可选 cfg、quality
% setMeshOptions:   payload.cfg；参数改变保留模型、清网格及下游
% setBoundary:      payload.cfg（整个最终配置）
% confirmBoundary:  可选 payload.freestream（按当前配置计算的来流）
% selectStep:       payload.step=1..5，只切换可用页面
% beginSolve:       无载荷；读取 state.frozenInputs 交给后端
% completeSolve:    payload.result、sourceRevision、solveId；可选 freestream
% failSolve:        payload.error、sourceRevision、solveId
%
% 所有 V 均已转换为米，F 为真实三角面索引。本函数不转换 STL 单位、
% 不猜测朝向、不添加压力/热流/温度字段，也不把输入壁温改成求解值。

if nargin<2,state=[];end
if nargin<3 || isempty(payload),payload=struct();end
assertScalarStruct(payload,'payload');
if ~(ischar(action) || (isstring(action) && isscalar(action)))
    error('RocketAeroHeat:WorkflowAction','工作流动作必须为单个文本。');
end
action=lower(char(action));

%% 模块一：建立尚未导入模型的状态，不默认为边界或结果已完成。
if strcmp(action,'init')
    cfg=optional(payload,'cfg',struct());assertScalarStruct(cfg,'cfg');
    state=struct('schemaVersion',1,'cfg',cfg,'rawMesh',[], ...
        'orientedMesh',[],'mesh',[],'quality',[],'modelSource',struct(), ...
        'freestream',[],'result',[],'researchResult',[],'revision',uint64(0), ...
        'modelRevision',uint64(0),'meshRevision',uint64(0), ...
        'boundaryRevision',uint64(0),'resultSourceRevision',NaN, ...
        'boundaryConfirmed',false,'currentStep',1,'busy',false, ...
        'status','awaiting-model','error',struct('identifier','','message',''), ...
        'solveId',uint64(0),'solveInputRevision',NaN,'frozenInputs',[], ...
        'lastInvalidation','初始化：尚未输入模型');
    state=deriveSteps(state);return
end
validateState(state);

switch action
    %% 模块二：模型或实际朝向改变，清除旧网格及下游结果。
    case 'invalidatemodel'
        cfg=optional(payload,'cfg',state.cfg);assertScalarStruct(cfg,'cfg');
        cfg=clearFaceBoundaryFields(cfg);
        state=invalidate(state,'文件或尺度已编辑，需重新输入并确认模型');
        state.revision=state.revision+1;state.modelRevision=state.modelRevision+1;
        state.cfg=cfg;state.rawMesh=[];state.orientedMesh=[];
        state.mesh=[];state.quality=[];state.modelSource=struct();
        state.currentStep=1;state.status='awaiting-model';
    case 'invalidateorientation'
        requireCondition(~isempty(state.rawMesh),'请先输入模型。');
        cfg=optional(payload,'cfg',state.cfg);assertScalarStruct(cfg,'cfg');
        cfg=clearFaceBoundaryFields(cfg);
        state=invalidate(state,'朝向预览已改变，需重新确认实际朝向');
        state.revision=state.revision+1;state.modelRevision=state.modelRevision+1;
        state.cfg=cfg;state.orientedMesh=[];state.mesh=[];state.quality=[];
        state.currentStep=1;state.status='model-ready';
    case 'setmodel'
        raw=required(payload,'rawMesh');validateMesh(raw,'rawMesh');
        raw=clearComputedFields(raw);
        oriented=optional(payload,'orientedMesh',[]);
        if ~isempty(oriented)
            validateMesh(oriented,'orientedMesh');oriented=clearComputedFields(oriented);
        end
        cfg=optional(payload,'cfg',state.cfg);assertScalarStruct(cfg,'cfg');
        source=optional(payload,'modelSource',struct());
        assertScalarStruct(source,'modelSource');
        if isfield(source,'units') && ~strcmp(char(source.units),'m')
            error('RocketAeroHeat:WorkflowUnits','工作流几何必须已换算为 SI 米。');
        end
        changed=~sameGeometry(raw,state.rawMesh) || ...
            ~sameGeometry(oriented,state.orientedMesh) || ...
            ~isequaln(cfg,state.cfg) || ~isequaln(source,state.modelSource);
        if changed
            state=invalidate(state,'模型、尺度或输入来源改变');
            cfg=clearFaceBoundaryFields(cfg);
            state.revision=state.revision+1;state.modelRevision=state.modelRevision+1;
            state.rawMesh=raw;state.orientedMesh=oriented;state.cfg=cfg;
            state.modelSource=source;state.mesh=[];state.quality=[];
            state.currentStep=1;state.status='model-ready';
        end
    case 'setorientation'
        requireCondition(~isempty(state.rawMesh),'请先输入模型。');
        oriented=required(payload,'orientedMesh');validateMesh(oriented,'orientedMesh');
        oriented=clearComputedFields(oriented);
        cfg=optional(payload,'cfg',state.cfg);assertScalarStruct(cfg,'cfg');
        changed=~sameGeometry(oriented,state.orientedMesh) || ~isequaln(cfg,state.cfg);
        if changed
            state=invalidate(state,'实际模型朝向或尺度改变');
            cfg=clearFaceBoundaryFields(cfg);
            state.revision=state.revision+1;state.modelRevision=state.modelRevision+1;
            state.orientedMesh=oriented;state.cfg=cfg;state.mesh=[];state.quality=[];
            state.currentStep=1;state.status='model-ready';
        end

    %% 模块三：确认真实网格；网格变化后保留配置，边界需重新审核。
    case 'setmeshoptions'
        cfg=required(payload,'cfg');assertScalarStruct(cfg,'cfg');
        requireUnchangedConfigSection(cfg,state.cfg,'geometry','模型输入或朝向');
        if ~isequaln(cfg,state.cfg)
            state=invalidate(state,'网格参数已改变，需重新生成网格并审核边界');
            cfg=clearFaceBoundaryFields(cfg);
            state.revision=state.revision+1;state.meshRevision=state.meshRevision+1;
            state.cfg=cfg;state.mesh=[];state.quality=[];
            state.currentStep=2;state.status='model-ready';
        end
    case 'setmesh'
        requireCondition(~isempty(state.orientedMesh),'请先确认模型朝向。');
        mesh=required(payload,'mesh');validateMesh(mesh,'mesh');
        mesh=clearComputedFields(mesh);
        cfg=optional(payload,'cfg',state.cfg);assertScalarStruct(cfg,'cfg');
        requireUnchangedConfigSection(cfg,state.cfg,'geometry','模型输入或朝向');
        quality=optional(payload,'quality',[]);
        changed=~isequaln(mesh,state.mesh) || ~isequaln(cfg,state.cfg);
        if changed
            state=invalidate(state,'求解网格或网格配置改变，边界条件需重新审核');
            cfg=clearFaceBoundaryFields(cfg);
            state.revision=state.revision+1;state.meshRevision=state.meshRevision+1;
            state.mesh=mesh;state.cfg=cfg;state.quality=quality;
            state.currentStep=2;state.status='mesh-ready';
        else
            state.quality=quality;
        end

    %% 模块四：应用边界配置和审核；给定壁温只存于输入配置。
    case 'setboundary'
        requireCondition(~isempty(state.mesh),'请先生成网格。');
        cfg=required(payload,'cfg');assertScalarStruct(cfg,'cfg');
        requireUnchangedConfigSection(cfg,state.cfg,'geometry','模型输入或朝向');
        requireUnchangedConfigSection(cfg,state.cfg,'remesh','网格参数');
        requireUnchangedConfigSection(cfg,state.cfg,'mesh','网格参数');
        if ~isequaln(cfg,state.cfg)
            state=invalidate(state,'边界条件或工程算法改变');
            state.revision=state.revision+1;state.boundaryRevision=state.boundaryRevision+1;
            state.cfg=cfg;state.currentStep=3;state.status='mesh-ready';
        end
    case 'confirmboundary'
        requireCondition(~isempty(state.mesh),'请先生成网格。');
        validateFlight(state.cfg);
        fs=optional(payload,'freestream',state.freestream);
        if ~isempty(fs),assertScalarStruct(fs,'freestream');end
        if ~isequaln(fs,state.freestream)
            state=invalidate(state,'确认当前边界及重新建立来流');
            state.revision=state.revision+1;state.boundaryRevision=state.boundaryRevision+1;
            state.freestream=fs;
        end
        state.boundaryConfirmed=true;
        if isempty(state.result) && isempty(state.researchResult) && ~state.busy,state.status='boundary-ready';end
        state.currentStep=3;

    %% 模块五：切换页面只读已有状态，不隐式启动计算。
    case 'selectstep'
        step=required(payload,'step');
        if ~isnumeric(step) || ~isscalar(step) || ~isfinite(step) || ...
                step~=fix(step) || step<1 || step>5
            error('RocketAeroHeat:WorkflowStep','工作流步骤必须为 1~5 的整数。');
        end
        state=deriveSteps(state);
        requireCondition(state.steps(step).available,'该步骤尚缺前置条件。');
        state.currentStep=double(step);

    %% 模块六：冻结求解输入，并对返回快照的版本、编号及几何严格守卫。
    case 'beginsolve'
        requireCondition(~state.busy,'已有当前求解正在运行。');
        requireCondition(~isempty(state.mesh) && state.boundaryConfirmed, ...
            '计算需要真实网格和已确认的边界条件。');
        validateFlight(state.cfg);
        state=clearResult(state);state.solveId=state.solveId+1;
        state.solveInputRevision=state.revision;
        state.frozenInputs=struct('mesh',state.mesh,'cfg',state.cfg, ...
            'freestream',state.freestream,'revision',state.revision,'solveId',state.solveId);
        state.busy=true;state.status='solving';state.currentStep=4;
        state.error=struct('identifier','','message','');
    case 'completesolve'
        validateSolveToken(state,payload);
        result=required(payload,'result');validateMesh(result,'result');
        if ~sameGeometry(result,state.mesh)
            error('RocketAeroHeat:WorkflowResultGeometry', ...
                '返回结果的 V/F 与当前求解网格不一致，拒绝接收。');
        end
        validateResultFields(result);
        % 不补场量或温度：后端快照中的 NaN、算法及逐面元状态完整保留。
        fs=optional(payload,'freestream',state.frozenInputs.freestream);
        if ~isempty(fs),assertScalarStruct(fs,'freestream');end
        state.result=result;state.freestream=fs;state.resultSourceRevision=state.revision;
        state.busy=false;state.status='solved';state.currentStep=5;
        state.frozenInputs=[];state.error=struct('identifier','','message','');
    case 'completeresearch'
        % 研究快照独立保存；不能通过普通completeSolve混入工程已计算面积。
        validateSolveToken(state,payload);
        result=required(payload,'result');validateMesh(result,'researchResult');
        requireCondition(sameGeometry(result,state.mesh),'研究场与当前原面元几何不一致。');
        requireCondition(isequal(result.productionQualified,false) && ...
            isequal(result.eligibleForSurfaceCoverage,false),'研究结果不能自授物理资格。');
        for name={'pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw','conservativeHeatLoadDensityRaw'}
            values=result.(name{1});
            requireCondition(isnumeric(values) && isvector(values) && numel(values)==size(state.mesh.F,1) && ...
                all(isfinite(values)),'研究场缺少完整原面元数值。');
        end
        requireCondition(strcmp(result.receipt.status,'terminal') && result.receipt.exitCode==0 && ...
            ~result.receipt.timedOut && ~result.receipt.cancelled,'研究快照需要实际成功进程终态。');
        state.researchResult=result;state.result=[];state.resultSourceRevision=state.revision;
        state.busy=false;state.status='research-ready';state.currentStep=5;state.frozenInputs=[];
        state.error=struct('identifier','','message','');
    case 'failsolve'
        validateSolveToken(state,payload);
        err=required(payload,'error');
        state=clearResult(state);state.busy=false;state.frozenInputs=[];
        state.status='failed';state.currentStep=4;
        state.error=errorRecord(err);
    otherwise
        error('RocketAeroHeat:WorkflowAction','未知工作流动作：%s。',action);
end
state=deriveSteps(state);
end

%% 模块七：集中派生侧边栏状态，结果完成不等于全面积物理模型已完成。
function state=deriveSteps(state)
modelReady=~isempty(state.orientedMesh);
meshReady=~isempty(state.mesh);
boundaryReady=meshReady && state.boundaryConfirmed;
resultReady=~isempty(state.result) && ...
    isequal(state.resultSourceRevision,state.revision) && ~state.busy;
researchReady=~isempty(state.researchResult) && ...
    isequal(state.resultSourceRevision,state.revision) && ~state.busy;
available=[true modelReady meshReady boundaryReady resultReady||researchReady];
completed=[modelReady meshReady boundaryReady resultReady resultReady];
keys={'model','mesh','boundary','solve','results'};
labels={'输入模型','生成网格','添加边界条件','计算','查看结果'};
template=struct('key','','label','','available',false,'completed',false,'status','pending');
state.steps=repmat(template,1,5);
for i=1:5
    state.steps(i).key=keys{i};state.steps(i).label=labels{i};
    state.steps(i).available=available(i);state.steps(i).completed=completed(i);
    if completed(i),state.steps(i).status='complete';
    elseif available(i),state.steps(i).status='available';end
end
if state.busy,state.steps(4).status='busy';end
if strcmp(state.status,'failed'),state.steps(4).status='failed';end
if researchReady,state.steps(4).status='research';state.steps(5).status='research';end
end

%% 模块八：所有输入失效都从同一入口清理，旧求解编号不复用。
function state=invalidate(state,reason)
hadFaceWall=isfield(state.cfg,'heating') && ...
    isfield(state.cfg.heating,'wallTemperatureByFace');
state=clearResult(state);state.freestream=[];state.boundaryConfirmed=false;
state.busy=false;state.frozenInputs=[];state.solveInputRevision=NaN;
state.error=struct('identifier','','message','');state.lastInvalidation=reason;
if hadFaceWall && ~strcmp(reason,'边界条件或工程算法改变') && ...
        ~strcmp(reason,'确认当前边界及重新建立来流')
    state.lastInvalidation=[reason '；旧逐面元壁温绑定已清除，需重新指定'];
end
end
function state=clearResult(state)
state.result=[];state.researchResult=[];state.resultSourceRevision=NaN;
end
function validateSolveToken(state,payload)
rev=required(payload,'sourceRevision');id=required(payload,'solveId');
if ~validCounter(rev) || ~validCounter(id) || ~state.busy || ...
        ~isequal(uint64(rev),state.revision) || ...
        ~isequal(uint64(rev),uint64(state.solveInputRevision)) || ...
        ~isequal(uint64(id),state.solveId)
    error('RocketAeroHeat:WorkflowStaleResult', ...
        '求解输入版本或求解编号已过期，拒绝将旧结果作为当前工况。');
end
end
function yes=validCounter(value)
yes=isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value>=0 && value==fix(value);
end

%% 模块九：验证真实 SI 三角网格；不依据尺寸猜单位或人为修补几何。
function validateMesh(mesh,label)
assertScalarStruct(mesh,label);
if ~isfield(mesh,'V') || ~isfield(mesh,'F')
    error('RocketAeroHeat:WorkflowMesh','%s 缺少真实 V/F。',label);
end
v=mesh.V;f=mesh.F;
if ~isnumeric(v) || ~isreal(v) || size(v,2)~=3 || size(v,1)<3 || ...
        ~ismatrix(v) || any(~isfinite(v(:))) || ...
        ~isnumeric(f) || ~isreal(f) || size(f,2)~=3 || isempty(f) || ...
        ~ismatrix(f) || any(~isfinite(f(:))) || any(f(:)~=fix(f(:))) || ...
        any(f(:)<1) || any(f(:)>size(v,1))
    error('RocketAeroHeat:WorkflowMesh','%s 必须含有限实数 V 和合法三角面 F。',label);
end
a=double(v(f(:,2),:))-double(v(f(:,1),:));
b=double(v(f(:,3),:))-double(v(f(:,1),:));
c=cross(a,b,2);
if any(all(c==0,2)) || any(~isfinite(c(:)))
    error('RocketAeroHeat:WorkflowMesh','%s 存在退化或非有限三角面。',label);
end
end
function yes=sameGeometry(a,b)
if isempty(a) || isempty(b),yes=isempty(a) && isempty(b);return;end
yes=isequaln(a.V,b.V) && isequaln(a.F,b.F);
end
function mesh=clearComputedFields(mesh)
% 输入网格可能来自先前快照；只保留几何、坐标、质量等上游信息。
names={'pressure','aerodynamicLoads','stagnation','boundaryLayerEdge', ...
    'boundaryLayer','heating','heatingCoverage','stagnationHeating', ...
    'wallThermal','aerothermalExposure','flightHeatingHistory', ...
    'surfaceStreamlines','axisymmetricMarch','externalFlow','baseFlow', ...
    'basePressure','baseCorrelation','sphericalNoseHeating'};
present=names(isfield(mesh,names));if ~isempty(present),mesh=rmfield(mesh,present);end
end
function cfg=clearFaceBoundaryFields(cfg)
% 逐面元输入随几何/编号绑定；同面数的新几何也不能继承旧向量。
if isfield(cfg,'heating') && isstruct(cfg.heating) && ...
        isfield(cfg.heating,'wallTemperatureByFace')
    cfg.heating=rmfield(cfg.heating,'wallTemperatureByFace');
end
end
function requireUnchangedConfigSection(a,b,key,label)
if ~isequaln(optional(a,key,[]),optional(b,key,[]))
    error('RocketAeroHeat:WorkflowInputHierarchy', ...
        '%s已改变，请通过对应上游步骤重新确认，不能与旧网格跨级提交。',label);
end
end
function validateResultFields(result)
% 仅真实计算快照可完成计算步骤；全 NaN 仍按未解状态原样保存。
count=size(result.F,1);groups={'pressure','pressure','pressure','heating'};
names={'staticPressure','Cp','pressureIncrement','heatFlux'};
for i=1:numel(names)
    group=groups{i};name=names{i};
    if ~isfield(result,group) || ~isstruct(result.(group)) || ...
            ~isscalar(result.(group)) || ~isfield(result.(group),name)
        error('RocketAeroHeat:WorkflowResultFields', ...
            '当前气动热快照缺少 %s.%s，不能标为已计算。',group,name);
    end
    value=result.(group).(name);
    if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
            numel(value)~=count || any(isinf(value(:)))
        error('RocketAeroHeat:WorkflowResultFields', ...
            '%s.%s 必须逐真实面元保存有限值或未解 NaN。',group,name);
    end
end
p=result.pressure;h=result.heating;
if ~isfield(p,'fieldStatus') || ~isnumeric(p.fieldStatus) || ...
        ~isreal(p.fieldStatus) || ~isvector(p.fieldStatus) || numel(p.fieldStatus)~=count || ...
        any(~ismember(p.fieldStatus(:),[0 1 2]))
    error('RocketAeroHeat:WorkflowResultFields', ...
        '压力必须逐面元记录 fieldStatus=0未解/1计算/2给定。');
end
missing=p.fieldStatus(:)==0;
for name={'staticPressure','Cp','pressureIncrement'}
    value=p.(name{1})(:);
    if any(isfinite(value(missing))) || any(~isfinite(value(~missing)))
        error('RocketAeroHeat:WorkflowResultFields', ...
            '压力 %s 的有限/NaN值与逐面元未解状态不一致。',name{1});
    end
end
if isfield(p,'missingPressureFace') && ...
        (~isequal(size(p.missingPressureFace),size(p.fieldStatus)) || ...
        ~isequal(logical(p.missingPressureFace(:)),missing))
    error('RocketAeroHeat:WorkflowResultFields','压力缺失掩码与 fieldStatus 不一致。');
end
if ~isfield(h,'eligibleFace') || ...
        ~(isnumeric(h.eligibleFace) || islogical(h.eligibleFace)) || ...
        ~isreal(h.eligibleFace) || ~isvector(h.eligibleFace) || numel(h.eligibleFace)~=count || ...
        any(~ismember(h.eligibleFace(:),[0 1])) || ...
        ~isequal(logical(h.eligibleFace(:)),isfinite(h.heatFlux(:)))
    error('RocketAeroHeat:WorkflowResultFields', ...
        '热流 eligibleFace 与真实热流有限值/未解NaN必须逐面元一致。');
end
end

%% 模块十：边界输入基础检查；完整物理适用性仍由后端各模型审核。
function validateFlight(cfg)
if ~isfield(cfg,'flight') || ~isstruct(cfg.flight) || ~isscalar(cfg.flight)
    error('RocketAeroHeat:WorkflowBoundary','尚未提供 cfg.flight 边界条件。');
end
names={'altitude','Mach','alpha','beta','Twall'};
for i=1:numel(names)
    name=names{i};
    if ~isfield(cfg.flight,name) || ~isnumeric(cfg.flight.(name)) || ...
            ~isreal(cfg.flight.(name)) || ~isscalar(cfg.flight.(name)) || ...
            ~isfinite(cfg.flight.(name))
        error('RocketAeroHeat:WorkflowBoundary','飞行输入 %s 必须为有限实数标量。',name);
    end
end
if cfg.flight.Mach<=0 || cfg.flight.Twall<=0
    error('RocketAeroHeat:WorkflowBoundary','Mach 和给定壁温必须大于零。');
end
end
function record=errorRecord(err)
if isa(err,'MException')
    record=struct('identifier',err.identifier,'message',err.message);
elseif isstruct(err) && isscalar(err) && isfield(err,'message')
    record=struct('identifier',optional(err,'identifier',''),'message',char(err.message));
elseif ischar(err) || (isstring(err) && isscalar(err))
    record=struct('identifier','','message',char(err));
else
    error('RocketAeroHeat:WorkflowError','失败必须保留异常或明确的错误文字。');
end
if isempty(strtrim(record.message))
    error('RocketAeroHeat:WorkflowError','不能使用空文字标记计算失败。');
end
end
function validateState(state)
assertScalarStruct(state,'state');
if ~isfield(state,'schemaVersion') || state.schemaVersion~=1
    error('RocketAeroHeat:WorkflowState','请先通过 init 建立工作流状态。');
end
end
function assertScalarStruct(value,label)
if ~isstruct(value) || ~isscalar(value)
    error('RocketAeroHeat:WorkflowPayload','%s 必须为标量结构。',label);
end
end
function value=required(payload,name)
if ~isfield(payload,name)
    error('RocketAeroHeat:WorkflowPayload','动作载荷缺少 %s。',name);
end
value=payload.(name);
end
function value=optional(payload,name,fallback)
if isfield(payload,name),value=payload.(name);else,value=fallback;end
end
function requireCondition(condition,message)
if ~condition,error('RocketAeroHeat:WorkflowPrerequisite','%s',message);end
end
