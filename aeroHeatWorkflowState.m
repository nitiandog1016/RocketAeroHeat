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
        if strcmp(optional(cfg,'computationBackend',''),'shortSupport65')
            [~,grid]=shortSupportImportBoundary(cfg.shortSupportBoundary.branch,optional(cfg,'shortSupportRadiusD',14));
            assert(isfield(mesh,'shortSupportImport') && isequal(mesh.shortSupportImport.source,state.modelSource) && ...
                strcmp(mesh.shortSupportImport.meshSHA256,grid.meshSHA256) && sameGeometry(mesh,state.orientedMesh), ...
                'RAH:Short:ImportMesh','请先核对当前 STL 原面和远场网格。');
        end
        if strcmp(optional(cfg,'computationBackend',''),'thermistorB') && isfield(state.modelSource,'file')
            assert(isfield(mesh,'thermistorImport') && isequal(mesh.thermistorImport.source,state.modelSource) && ...
                sameGeometry(mesh,state.orientedMesh),'RAH:Bead:ImportMesh','请用当前 STL 文件和确认朝向生成试件体网。');
        end
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
    case 'invalidateboundary'
        % 界面正在编辑但验证失败时也须撤销旧确认，防止以旧边界悄悄运行。
        requireCondition(~state.busy && ~isempty(state.mesh),'只能使当前空闲网格的边界失效。');
        state=invalidate(state,'边界编辑未通过验证，须修正并重新确认');
        state.revision=state.revision+1;state.boundaryRevision=state.boundaryRevision+1;
        state.currentStep=3;state.status='mesh-ready';
    case 'setboundary'
        requireCondition(~isempty(state.mesh),'请先生成网格。');
        cfg=required(payload,'cfg');assertScalarStruct(cfg,'cfg');
        requireUnchangedConfigSection(cfg,state.cfg,'geometry','模型输入或朝向');
        requireUnchangedConfigSection(cfg,state.cfg,'remesh','网格参数');
        requireUnchangedConfigSection(cfg,state.cfg,'mesh','网格参数');
        if any(strcmp({optional(state.cfg,'computationBackend',''),optional(cfg,'computationBackend','')},'shortSupport65'))
            requireUnchangedConfigSection(cfg,state.cfg,'computationBackend','银球体网');
            requireUnchangedConfigSection(cfg,state.cfg,'shortSupportRadiusD','银球远场半径');
        end
        if strcmp(optional(state.cfg,'computationBackend',''),'thermistorB') || strcmp(optional(cfg,'computationBackend',''),'thermistorB')
            requireUnchangedConfigSection(cfg,state.cfg,'computationBackend','试件来源');
            requireUnchangedConfigSection(cfg,state.cfg,'thermistorAngularCells','周向网格');
            requireUnchangedConfigSection(cfg,state.cfg,'thermistorVariant','试件外形');
            requireUnchangedConfigSection(cfg,state.cfg,'thermistorRefineLeads','带线中段网距');
        end
        if ~isequaln(cfg,state.cfg)
            state=invalidate(state,'边界条件或工程算法改变');
            state.revision=state.revision+1;state.boundaryRevision=state.boundaryRevision+1;
            state.cfg=cfg;state.currentStep=3;state.status='mesh-ready';
        end
    case 'confirmboundary'
        requireCondition(~isempty(state.mesh),'请先生成网格。');
        validateFlight(state.cfg);
        fs=optional(payload,'freestream',state.freestream);
        if strcmp(optional(state.cfg,'computationBackend',''),'shortSupport65')
            fs=optional(payload,'freestream',state.cfg.shortSupportBoundary);
            assert(isequal(fs,state.cfg.shortSupportBoundary),'RAH:Short:Boundary','来流须对应所确认的温界和步长。');
        end
        if strcmp(optional(state.cfg,'computationBackend',''),'thermistorB')
            [~,reference]=thermistorWorkbenchBoundary(state.cfg.thermistorBoundary,state.cfg.thermistorAngularCells,state.cfg);
            if isempty(fs),fs=reference;end
            assert(isequal(fs,reference),'RAH:Bead:Boundary','来流记录与确认的珠体边界不符。');
        end
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
        if strcmp(optional(state.cfg,'computationBackend',''),'shortSupport65')
            % 动作与目录在按下按钮时冻结。仅新算生成独立任务号，读档无需本机体场。
            source=validatestring(required(payload,'shortSupportAction'),{'archive','solve'});
            folder='';
            if strcmp(source,'solve')
                folder=fullfile(fileparts(mfilename('fullpath')),'outputs',['short65_time_' char(java.util.UUID.randomUUID())]);
            end
            state.frozenInputs.cfg.shortSupportAction=source;state.frozenInputs.cfg.shortSupportJob=folder;
        elseif strcmp(optional(state.cfg,'computationBackend',''),'thermistorB') && state.cfg.thermistorAngularCells==48 && ~strcmp(optional(state.cfg,'thermistorVariant','smooth'),'sidewind')
            % 将动作与目录一起冻结，避免同一几何和同刻的其他计算被误接收。
            source=validatestring(optional(payload,'startupAction','run'),{'run','read','resume'});
            if strcmp(source,'run')
                folder=fullfile(fileparts(mfilename('fullpath')),'outputs','v294_validation','workbench',char(java.util.UUID.randomUUID()));
            else
                folder=required(payload,'runDirectory');
                requireCondition(isfolder(folder),'请选择存在的计算根目录。');
                folder=char(java.io.File(folder).getCanonicalPath());
            end
            state.frozenInputs.cfg.thermistorResultSource=source;
            state.frozenInputs.cfg.thermistorRunDirectory=folder;
        elseif strcmp(optional(state.cfg,'computationBackend',''),'thermistorB') && any(strcmp(optional(state.cfg,'thermistorVariant','smooth'),{'leads','midbody','sidewind'}))
            % 在按下计算按钮时确定同一任务的来源，防止另一次同工况的结果串入。
            measured=strcmp(state.cfg.thermistorVariant,'sidewind') && state.cfg.thermistorAngularCells==96;
            wire=measured && isfield(state.cfg.thermistorBoundary,'leadThermalModel');
            reading=wire && isfield(state.cfg.thermistorBoundary,'resistanceTemperatureK');
            choices={'solve','archive'};
            if wire && ~reading,choices=[choices {'resumeWire'}];elseif measured && ~wire,choices=[choices {'cold','resume'}];end
            source=validatestring(optional(payload,'leadAction','solve'),choices);folder='';
            if wire && ~reading && ~strcmp(source,'archive')
                state.frozenInputs.cfg.thermistorInvocation=char(java.util.UUID.randomUUID());
                if strcmp(source,'resumeWire')
                    folder=required(payload,'runDirectory');
                    requireCondition(isfile(fullfile(folder,'wire_job.json')) || isfile(fullfile(folder,'prepared','profile.json')) && ...
                        isfile(fullfile(folder,'completed.json')),'请选择长引线任务目录或已完成的原生计算目录。');
                    folder=char(java.io.File(folder).getCanonicalPath());
                else
                    folder=fullfile(fileparts(mfilename('fullpath')),'outputs',['thermistor_wire_' char(java.util.UUID.randomUUID())]);
                    seed=optional(payload,'wireSeedDirectory','');
                    requireCondition(isempty(seed) || isfile(fullfile(seed,'completed.json')),'均温初猜请选择已完成的 96 档冷启动任务。');
                    state.frozenInputs.cfg.wireSeedDirectory=seed;
                end
            elseif measured && any(strcmp(source,{'cold','resume'}))
                token=char(java.util.UUID.randomUUID());
                state.frozenInputs.cfg.thermistorInvocation=token;
                if strcmp(source,'resume')
                    folder=required(payload,'runDirectory');
                    requireCondition(isfile(fullfile(folder,'job.json')),'请选择有冷启动任务记录的根目录。');
                    folder=char(java.io.File(folder).getCanonicalPath());
                else
                    folder=fullfile(fileparts(mfilename('fullpath')),'outputs',['thermistor_measured_' token]);
                end
            elseif strcmp(source,'solve') && ~measured
                folder=fullfile(fileparts(mfilename('fullpath')),'outputs',['thermistor_' state.cfg.thermistorVariant '_run_' char(java.util.UUID.randomUUID())]);
            end
            state.frozenInputs.cfg.thermistorResultSource=source;
            state.frozenInputs.cfg.thermistorRunDirectory=folder;
        end
        state.busy=true;state.status='solving';state.currentStep=4;
        state.error=struct('identifier','','message','');
    case 'completesolve'
        validateSolveToken(state,payload);
        result=required(payload,'result');validateMesh(result,'result');
        if ~sameGeometry(result,state.mesh)
            error('RocketAeroHeat:WorkflowResultGeometry', ...
                '返回结果的 V/F 与当前求解网格不一致，拒绝接收。');
        end
        if strcmp(optional(state.cfg,'computationBackend',''),'thermistorB')
            validateThermistorResult(result,state.frozenInputs);
        elseif strcmp(optional(state.cfg,'computationBackend',''),'shortSupport65')
            cfg=state.frozenInputs.cfg;b=cfg.shortSupportBoundary;binding=state.mesh.shortSupportImport;
            w=result.hybridWall;id=result.evidenceIdentity;n=size(result.F,1);
            assert(isequal(result.shortSupportImport,binding) && strcmp(id.sourceFrame.stlSHA256,binding.source.sha256) && ...
                strcmp(id.meshSHA256,binding.meshSHA256) && strcmp(id.sourceFrame.caseSHA256,b.caseSHA256) && ...
                strcmp(id.sourceFrame.initialSHA256,b.initialFrame) && strcmp(result.shortSupportAction,cfg.shortSupportAction) && ...
                strcmp(result.shortSupportJob,cfg.shortSupportJob) && strcmp(result.shortSupportDirectionGrid,['time' b.branch]) && ...
                result.hybridCondition.elapsedPhysicalTimeS==b.elapsedS && result.hybridCondition.physicalStepSizeS==b.dtS && ...
                isequal(required(payload,'freestream'),b) && isequal(w.faceID,(0:n-1)') && ...
                all(isfinite(w.pressurePa) & isfinite(w.heatFluxWm2) & isfinite(w.nearWallCellTemperatureK)) && ...
                numel(result.hybridWallTemperatureK)==n && all(isfinite(result.hybridWallTemperatureK)), ...
                'RAH:Short:Snapshot','返回的银球快照与本次冻结输入或原面场不符。');
            if strcmp(cfg.shortSupportAction,'archive')
                [~,pair]=shortSupportImportBoundary(b.branch,optional(cfg,'shortSupportRadiusD',14));
                assert(strcmp(id.sourceFrame.sha256,pair.([b.branch 'Frame'])),'RAH:Short:Snapshot','归档帧与所选步长不符。');
            end
        elseif strcmp(optional(state.cfg,'computationBackend',''),'cm011Imported')
            validateCM011Result(result,state.mesh);
        else
            validateResultFields(result);
        end
        % 不补场量或温度：后端快照中的 NaN、算法及逐面元状态完整保留。
        fs=optional(payload,'freestream',state.frozenInputs.freestream);
        if strcmp(optional(state.cfg,'computationBackend',''),'thermistorB')
            assert(isequal(fs,state.frozenInputs.freestream),'RAH:Bead:Snapshot','返回的来流与本次冻结工况不符。');
        end
        if ~isempty(fs),assertScalarStruct(fs,'freestream');end
        state.result=result;state.freestream=fs;state.resultSourceRevision=state.revision;
        state.busy=false;state.status='solved';state.currentStep=5;
        state.frozenInputs=[];state.error=struct('identifier','','message','');
    case {'completeresearch','completecheckpointresearch'}
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
        if strcmp(action,'completecheckpointresearch')
            % 仅由显式超时重开入口接受；普通成功任务路径仍保持原有严格条件。
            requireCondition(strcmp(result.receipt.status,'terminal') && result.receipt.timedOut && ...
                ~result.receipt.cancelled && isfield(result,'explicitTimeoutReview') && ...
                isfield(result,'convergenceAudit') && result.convergenceAudit.checkpointIdentityConfirmed, ...
                '超时研究查看需要真实终态和已确认的原检查点身份。');
        else
            requireCondition(strcmp(result.receipt.status,'terminal') && result.receipt.exitCode==0 && ...
                ~result.receipt.timedOut && ~result.receipt.cancelled,'研究快照需要实际成功进程终态。');
        end
        state.researchResult=result;state.result=[];state.resultSourceRevision=state.revision;
        state.busy=false;state.status='research-ready';state.currentStep=5;state.frozenInputs=[];
        state.error=struct('identifier','','message','');
    case 'pausesolve'
        % 本次热敏珠原生操作已在保存点退出，边界仍可用于接续。
        validateSolveToken(state,payload);
        requireCondition(isfield(state.frozenInputs.cfg,'thermistorInvocation'),'当前任务没有冷启动保存点。');
        state=clearResult(state);state.busy=false;state.frozenInputs=[];
        state.status='paused';state.currentStep=4;
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
if strcmp(state.status,'paused'),state.steps(4).status='paused';end
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
if strcmp(optional(cfg,'computationBackend',''),'shortSupport65')
    b=cfg.shortSupportBoundary;
    assert(isequal(b,shortSupportImportBoundary(b.branch,optional(cfg,'shortSupportRadiusD',14))),'RAH:Short:Boundary','请使用当前双支管试件登记的工况和物理时间。');return;
end
if strcmp(optional(cfg,'computationBackend',''),'thermistorB')
    thermistorWorkbenchBoundary(cfg.thermistorBoundary,cfg.thermistorAngularCells,cfg);
    return;
end
if strcmp(optional(cfg,'computationBackend',''),'cm011Imported')
    cm011WorkbenchBoundary(cfg.cm011Boundary);return;
end
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
function validateThermistorResult(result,frozen)
% 核对低速试件的冻结边界，分别检查给定壁温和由净功率求解的稳态壁温。
cfg=frozen.cfg;[b,fs]=thermistorWorkbenchBoundary(cfg.thermistorBoundary,cfg.thermistorAngularCells,cfg);
p=frozen.mesh.thermistorPrepared;w=result.hybridWall;n=size(frozen.mesh.F,1);c=result.hybridCondition;
id=result.evidenceIdentity;
if isfield(frozen.mesh,'thermistorImport')
    map=frozen.mesh.thermistorImport;
    assert(isequal(map,result.thermistorImport) && strcmp(id.sourceFrame.stlSHA256,map.source.sha256) && ...
        strcmp(cm011FileSHA256(map.source.file),map.source.sha256),'RAH:Bead:Snapshot','STL 来源或逐面对应与当前工况不同。');
end
wire=isfield(b,'leadThermalModel');
if wire
    % 逐面给定的长引线热边界与本次动作、壁温分布指纹和全表面气动热源同时绑定。
    m=result.thermistorWire;
    assert(strcmp(m.action,cfg.thermistorResultSource) && strcmpi(id.runDirectory,cfg.thermistorRunDirectory) && ...
        id.reusedArchive==(isfield(b,'resistanceTemperatureK') || strcmp(m.action,'archive')) && strcmp(m.sourceFrameSHA256,id.sourceFrame.sha256) && ...
        strcmp(m.wallSHA256,id.sourceFrame.wallSHA256) && strcmp(c.boundarySHA256,m.boundarySHA256) && ...
        strcmp(c.leadThermalModel,b.leadThermalModel) && m.faces==n && m.jointExperimentVerifiedM2==0);
elseif strcmp(optional(cfg,'thermistorVariant','smooth'),'sidewind') && p.angularCells==96
    % 12 K 测点绑定当前动作、原面与源帧；冷启动还须绑定冻结目录及本次操作的回执。
    m=result.thermistorMeasuredPoint;src=m.scaledFrom;
    if strcmp(cfg.thermistorResultSource,'archive'),src=m.nativeCheck;end
    fresh=any(strcmp(cfg.thermistorResultSource,{'cold','resume'}));
    if fresh
        src=m.coldSource;
        assert(strcmp(src.invocation,cfg.thermistorInvocation) && strcmp(src.initialKind,'geometry-and-boundary'));
    end
    assert(strcmp(id.sourceFrame.sha256,src.frame) && strcmp(id.sourceFrame.wallSHA256,src.wallSHA256) && ...
        strcmp(id.meshSHA256,m.meshSHA256) && m.deltaTK==12 && m.wallK==b.wallTemperatureK && ...
        strcmp(m.action,cfg.thermistorResultSource) && strcmpi(id.runDirectory,cfg.thermistorRunDirectory) && id.reusedArchive==~fresh && ...
        abs(m.wholeW-sum([m.regions.outwardW]))<1e-12 && m.jointVerifiedM2==0,'RAH:Bead:Point','测点结果与边界或来源不符。');
elseif strcmp(optional(cfg,'thermistorVariant','smooth'),'sidewind')
    s=result.thermistorSidewind;
    assert(p.angularCells==48 && strcmp(s.meshSHA256,p.meshSHA256) && strcmp(s.sourceFrameSHA256,id.sourceFrame.sha256) && ...
        strcmp(s.wallSHA256,id.sourceFrame.wallSHA256) && numel(s.regions)==7 && sum([s.regions.faces])==n && ...
        strcmpi(id.runDirectory,cfg.thermistorRunDirectory) && id.reusedArchive==strcmp(cfg.thermistorResultSource,'archive'), ...
        'RAH:Lead:Snapshot','侧吹源帧、全周体网或区域与冻结工况不符。');
elseif p.angularCells==48
    % 全周计算核对两个实际物理半步及给定壁温；扇区稳态沿用原有资格。
    s=result.beadNearWallComparison.summary;
    assert(c.elapsedPhysicalTimeS==7.80936086544745e-6 && strcmp(c.initialStateSHA256,p.seedSHA256) && ...
        strcmp(s.fineSurface.meshSHA256,p.meshSHA256) && strcmp(s.fineSurface.energy.acceptedFrameSHA256,id.sourceFrame.sha256) && ...
        strcmp(s.firstHalfEnergy.acceptedFrameSHA256,s.fineSurface.energy.oldFrameSHA256) && ...
        strcmpi(id.runDirectory,cfg.thermistorRunDirectory) && strcmp(id.startupAction,cfg.thermistorResultSource), ...
        'RAH:Bead:Snapshot','计算目录、终点或源帧不符。');
elseif any(strcmp(optional(cfg,'thermistorVariant','smooth'),{'leads','midbody'}))
    % 按试件分别核对带线七区或光滑三区，冻结快照必须属于本次体网。
    if strcmp(cfg.thermistorVariant,'leads'),s=result.thermistorLead;groups=7;
    else,s=result.thermistorMidbody;groups=3;end
    native=s;if isfield(s,'referenceNative'),native=s.referenceNative;end
    assert(strcmp(native.meshSHA256,p.meshSHA256) && strcmp(id.sourceFrame.sha256,native.sourceFrameSHA256) && ...
        strcmp(id.sourceFrame.wallSHA256,native.wallSHA256) && numel(s.regions)==groups && sum([s.regions.faces])==n && ...
        strcmpi(id.runDirectory,cfg.thermistorRunDirectory) && id.reusedArchive==strcmp(cfg.thermistorResultSource,'archive'), ...
        'RAH:Lead:Snapshot','带线原面、区域或终场不符。');
else
    assert(result.thermistorBead.cellAngleDeg==30/p.angularCells,'RAH:Bead:Snapshot','扇区分辨率不符。');
end
assert(cfg.thermistorAngularCells==p.angularCells && ...
    strcmp(id.origin,'research') && isequal(id.productionQualified,false) && strcmp(id.meshSHA256,p.meshSHA256) && ...
    strcmp(id.preparedSHA256,p.manifestSHA256) && strcmp(id.sourceFrame.bodySHA256,p.bodySHA256), ...
    'RAH:Bead:Snapshot','热敏珠快照的网格或来源不符。');
assert(isequal(c.velocityMps,fs.velocityMps) && c.ambientPressurePa==b.airPressurePa && ...
    c.ambientTemperatureK==b.airTemperatureK && isequal(c,id.physicalCondition), ...
    'RAH:Bead:Snapshot','结果边界与当前确认的试件工况不符。');
if wire
    if isfield(b,'resistanceTemperatureK')
        t=m.thermometry;beadK=t.surfaceK;
        assert(c.resistanceTemperatureK==b.resistanceTemperatureK && t.internalK==c.resistanceTemperatureK && ...
            abs((beadK-c.ambientTemperatureK)*(1+t.ratio)-(t.internalK-c.ambientTemperatureK))<1e-10);
    else,beadK=b.wallTemperatureK;end
    tw=result.hybridWallTemperatureK;assert(~isfield(result,'hybridSolvedWallTemperatureK') && numel(tw)==n && ...
        c.beadWallK==beadK && abs(min(tw)-m.profile.endK)<1e-9 && abs(max(tw)-c.beadWallK)<1e-9);
elseif isfield(b,'netHeatingPowerW')
    if isfield(result,'thermistorLead'),t=result.thermistorLead.thermal;
    else,t=result.thermistorBead.thermal;end
    assert(isfield(c,'netHeatingPowerW') && c.netHeatingPowerW==b.netHeatingPowerW && ...
        ~isfield(c,'wallTemperatureK') && ~isfield(result,'hybridWallTemperatureK') && ...
        t.netHeatingPowerW==b.netHeatingPowerW && t.conductanceWK>0 && ...
        abs(t.wallTemperatureK-b.airTemperatureK-b.netHeatingPowerW/t.conductanceWK)<1e-10 && ...
        isequal(result.hybridSolvedWallTemperatureK,repmat(t.wallTemperatureK,n,1)) && ...
        abs(t.surfaceReleaseW-b.netHeatingPowerW)<=1e-11*b.netHeatingPowerW, ...
        'RAH:Bead:Snapshot','稳态壁温或空气侧放热与冻结净功率不符。');
else
    assert(~isfield(result,'hybridSolvedWallTemperatureK') && c.wallTemperatureK==b.wallTemperatureK && ...
        isequal(result.hybridWallTemperatureK,repmat(b.wallTemperatureK,n,1)), ...
        'RAH:Bead:Snapshot','给定壁温与冻结输入不符。');
end
assert(istable(w) && isequal(w.faceID,(0:n-1)') && all(isfinite(w.pressurePa) & w.pressurePa>0) && ...
    all(isfinite(w.heatFluxWm2)) && all(isfinite(w.nearWallCellTemperatureK)) && ...
    all(strlength(w.pressureMethod)>0 & strlength(w.heatFluxMethod)>0), ...
    'RAH:Bead:Snapshot','原面数值或算法来源不完整。');
end
function validateCM011Result(result,mesh)
% 低速全壁数值快照须与本次网格指纹及原面编号一致，各场量使用真实来源。
p=mesh.cm011Prepared;n=size(mesh.F,1);w=result.hybridWall;
assert(istable(w) && height(w)==n && isequal(w.faceID,(0:n-1)') && ...
    strcmp(result.evidenceIdentity.stlSHA256,p.sourceSHA256) && ...
    strcmp(result.evidenceIdentity.preparedSHA256,p.manifestSHA256) && ...
    isequal(result.evidenceIdentity.productionQualified,false),'RAH:CM011:Snapshot','试件快照的来源或面编号不符。');
for name={'pressurePa','pressureGaugePa','heatFluxWm2','wallFluidTemperatureK','nearWallCellTemperatureK'}
    x=w.(name{1});assert(isequal(size(x),[n,1]) && isreal(x) && all(isfinite(x)),'RAH:CM011:Snapshot','试件场量不完整。');
end
assert(all(w.pressurePa>0) && all(result.hybridWallTemperatureK==300.1) && ...
    all(strlength(w.pressureMethod)>0) && all(strlength(w.heatFluxMethod)>0),'RAH:CM011:Snapshot','压力、壁温或算法记录无效。');
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
