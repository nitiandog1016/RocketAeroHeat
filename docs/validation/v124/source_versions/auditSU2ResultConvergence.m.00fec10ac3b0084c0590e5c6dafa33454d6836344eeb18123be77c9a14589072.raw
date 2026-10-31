function audit=auditSU2ResultConvergence(folder,result)
%AUDITSU2RESULTCONVERGENCE 区分进程终态、历史末行和实际显示检查点，不推定物理资格。
% 模块一：读取原配置门槛和完整连续历史，按原进程回执核对文件SHA。
% 模块二：按原生面积热负荷唯一关联保存时刻；重启动与壁面p/T再做EOS核对。
% 模块三：缺证据、重复候选、未达残差门槛分别保留；不把退出0当作收敛。
audit=struct('schemaVersion',1,'status','evidence-unavailable','reason','缺少原配置或完整历史，不能推断收敛。', ...
    'terminalIteration',NaN,'checkpointIteration',NaN,'checkpointIdentityConfirmed',false, ...
    'terminalResidualPassed',false,'checkpointResidualPassed',false,'productionQualified',false, ...
    'historyRows',0,'threshold',NaN,'fields',{{}},'terminalResiduals',[],'checkpointResiduals',[], ...
    'last250HeatPowerRelativeChange',NaN,'evidenceSHA256',struct,'folder',folder);
try
    needed={'case.cfg','history.csv','surface_flow.vtk','engine_process.json'};
    if ~all(cellfun(@(p)isfile(fullfile(folder,p)),needed)),return;end
    receiptText=fileread(fullfile(folder,'engine_process.json'));receipt=jsondecode(receiptText);
    assert(strcmp(receipt.status,'terminal') && isequaln(receipt,result.receipt),'结果与真实终态回执不一致。');
    % MATLAB会截断过长JSON键，绝对路径可能碰撞；散列/稳定性须检查原始JSON键。
    stable=regexp(receiptText,'(?s)"inputsStable"\s*:\s*\{([^}]+)\}','tokens','once');
    assert(~isempty(stable) && isempty(regexp(stable{1},':\s*false','once')) && contains(stable{1},'true'),'运行时输入有变更。');
    for name={'history.csv','surface_flow.vtk'}
        key=matlab.lang.makeValidName(name{1});h=fileSHA(fullfile(folder,name{1}));
        assert(isfield(receipt.outputSHA256,key) && strcmp(h,receipt.outputSHA256.(key)),'输出SHA与进程记录不一致。');
        audit.evidenceSHA256.(key)=h;
    end
    config=fileread(fullfile(folder,'case.cfg'));configPath=fullfile(receipt.request.workDirectory,receipt.request.configFile);
    h=fileSHA(fullfile(folder,'case.cfg'));
    token=regexp(receiptText,[regexptranslate('escape',jsonencode(configPath)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
    assert(numel(token)==1 && strcmp(h,token{1}{1}),'配置SHA不匹配。');
    audit.evidenceSHA256.case_cfg=h;audit.evidenceSHA256.engine_process_json=fileSHA(fullfile(folder,'engine_process.json'));
    assert(isempty(regexp(config,'(?m)^\s*TIME_DOMAIN\s*=\s*YES','once')),'当前仅审核稳态历史。');
    threshold=number(config,'CONV_RESIDUAL_MINVAL');start=number(config,'CONV_STARTITER');
    cadence=number(config,'OUTPUT_WRT_FREQ');assert(cadence>=1 && cadence==fix(cadence));
    names=regexp(config,'(?m)^\s*CONV_FIELD\s*=\s*\(([^)]*)\)','tokens','once');assert(~isempty(names));
    selected=strtrim(strsplit(names{1},','));
    allowed={'RMS_DENSITY','RMS_MOMENTUM-X','RMS_MOMENTUM-Y','RMS_MOMENTUM-Z','RMS_ENERGY','RMS_TKE','RMS_DISSIPATION'};
    columns={'rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]','rms[k]','rms[w]'};
    assert(all(ismember(selected,allowed)),'包含尚未支持的非残差收敛条件。');
    fid=fopen(fullfile(folder,'history.csv'),'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));header=fgetl(fid);clear guard;
    keys=regexp(header,'"([^"]+)"','tokens');keys=cellfun(@(c)strtrim(c{1}),keys,'UniformOutput',false);
    history=readmatrix(fullfile(folder,'history.csv'),'NumHeaderLines',1);assert(size(history,2)==numel(keys));
    assert(~isempty(history),'历史为空。');
    % 固定引擎的二维Buffet诊断可能未定义；它不是收敛量或热负荷，明确单列。
    % 只允许该诊断的NaN，其他列（尤其残差/HF）及任意Inf仍必须拒绝。
    unused=strcmp(keys,'Buffet');requiredHistory=history(:,~unused);
    assert(all(isfinite(requiredHistory(:))) && ~any(isinf(history(:))),'历史含不完整/非有限必要字段。');
    audit.unavailableDiagnosticFields=keys(any(isnan(history),1));
    iteration=history(:,strcmp(keys,'Inner_Iter'));assert(~isempty(iteration) && iteration(1)==0 && all(diff(iteration)==1),'必须是完整连续历史，不能使用稀疏节选。');
    residual=zeros(size(history,1),numel(selected));
    for j=1:numel(selected),residual(:,j)=history(:,strcmp(keys,columns{strcmp(allowed,selected{j})}));end
    audit.fields=selected;audit.threshold=threshold;audit.historyRows=size(history,1);
    audit.terminalIteration=iteration(end);audit.terminalResiduals=residual(end,:);
    audit.terminalResidualPassed=iteration(end)>=start && all(residual(end,:)<threshold);
    hf=history(:,strcmp(keys,'HF'));assert(numel(hf)==numel(iteration));
    if numel(hf)>=251,audit.last250HeatPowerRelativeChange=(hf(end)-hf(end-250))/max(abs(hf(end-250)),1);end
    candidates=find(mod(iteration,cadence)==0);
    if receipt.exitCode==0 && ~receipt.timedOut && ~receipt.cancelled,candidates=unique([candidates;numel(iteration)]);end
    relative=abs(hf(candidates)-result.reconstructedEngineHeatPowerRaw)/max(abs(result.reconstructedEngineHeatPowerRaw),1);
    matches=candidates(relative<2e-6);audit.heatLoadCandidateCount=numel(matches);
    audit.status='history-audited-checkpoint-unconfirmed';audit.reason='保存负荷无法唯一对应输出时刻；显示迭代保持未知。';
    native=identifySU2NativeCheckpoint(folder,receipt,config);audit.nativeCheckpointIdentity=native;
    if native.enabled
        if ~native.confirmed,audit.reason=native.reason;return;end
        k=find(iteration==native.localIteration);assert(isscalar(k),'具名迭代不在完整历史中。');
        assert(abs(hf(k)-result.reconstructedEngineHeatPowerRaw)/max(abs(hf(k)),1)<2e-6,'具名场与同迭代热负荷不一致。');
        audit.checkpointIdentityMethod='native-indexed-pair-SHA-and-EOS';
    else
        if numel(matches)~=1,return;end
        k=matches(1);audit.checkpointIdentityMethod='unique-native-heat-load-and-EOS';
    end
    audit.checkpointIteration=iteration(k);audit.checkpointResiduals=residual(k,:);
    audit.checkpointResidualPassed=iteration(k)>=start && all(residual(k,:)<threshold);
    audit.heatLoadMatchRelativeError=abs(hf(k)-result.reconstructedEngineHeatPowerRaw)/max(abs(hf(k)),1);
    % 热负荷相近不够：保存的体场和壁面文件须同态；二维/三维均用实际列数处理。
    restart=fullfile(folder,'restart_flow.csv');
    if ~isfile(restart),audit.reason='负荷唯一关联，但缺少重启动文件，未确认场身份。';return;end
    h=fileSHA(restart);assert(strcmp(h,receipt.outputSHA256.restart_flow_csv),'重启动SHA不匹配。');
    audit.evidenceSHA256.restart_flow_csv=h;
    fid=fopen(restart,'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));header=fgetl(fid);clear guard;
    cols=regexp(header,'"([^"]+)"','tokens');cols=cellfun(@(x)x{1},cols,'UniformOutput',false);
    nd=2+any(strcmp(cols,'z'));expected={'PointID','x','y'};if nd==3,expected{end+1}='z';end
    expected=[expected,{'Density','Momentum_x','Momentum_y'}];if nd==3,expected{end+1}='Momentum_z';end
    expected{end+1}='Energy';sst=~isempty(regexp(config,'(?m)^\s*SOLVER\s*=\s*RANS\s*$','once'));
    if sst
        assert(isfield(receipt,'physicsModel') && strcmp(receipt.physicsModel.model,'SST-2003m'),'RANS缺少实际模型合同。');
        assert(all(ismember({'RMS_DENSITY','RMS_ENERGY','RMS_TKE','RMS_DISSIPATION'},selected)));
        expected=[expected,{'Turb_Kin_Energy','Omega'}];
    end
    assert(isequal(cols,expected),'重启动列名与模型不一致。');
    rest=readmatrix(restart,'NumHeaderLines',1);assert(size(rest,2)==numel(cols));
    assert(all(isfinite(rest(:))) && isequal(rest(:,1),(0:size(rest,1)-1)'),'重启动不完整。');
    rho=rest(:,nd+2);velocity=rest(:,nd+3:2*nd+2)./rho;internal=rest(:,strcmp(cols,'Energy'))./rho-sum(velocity.^2,2)/2;
    if sst
        k=rest(:,strcmp(cols,'Turb_Kin_Energy'));w=rest(:,strcmp(cols,'Omega'));
        assert(all(k>=0) && all(w>0),'湍流状态不满足k非负/ω正。');internal=internal-k;
        audit.physicsModel=receipt.physicsModel;audit.turbulentKRange=[min(k) max(k)];audit.omegaRange=[min(w) max(w)];
    end
    assert(all(rho>0));
    if contains(config,'FLUID_MODEL= IDEAL_GAS')
        assert(all(internal>0));gamma=number(config,'GAMMA_VALUE');R=number(config,'GAS_CONSTANT');
        temperature=internal*(gamma-1)/R;pressure=rho.*temperature*R;
    elseif contains(config,'FLUID_MODEL= COOLPROP')
        [pressure,temperature,thermo]=readSU2CoolPropStateAudit(folder,receipt,rest);
        audit.thermodynamicModel=thermo.thermodynamicModel;audit.thermodynamicAudit=thermo;
    else,error('当前EOS未实现独立状态审核，不能套用定比热公式。');end
    wall=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));points=[rest(:,2:nd+1) zeros(size(rest,1),3-nd)];
    ids=zeros(size(wall.vertices,1),1);distance=ids;
    for j=1:numel(ids),[d,ids(j)]=min(sum((points-wall.vertices(j,:)).^2,2));distance(j)=sqrt(d);end
    assert(numel(unique(ids))==numel(ids) && max(distance)<8e-7*max(1,max(abs(points(:)))),'壁节点匹配不唯一。');
    audit.wallEOSPressureError=max(abs(wall.pressureRaw./pressure(ids)-1));
    audit.wallEOSTemperatureError=max(abs(wall.temperatureRaw./temperature(ids)-1));
    assert(max([audit.wallEOSPressureError audit.wallEOSTemperatureError])<2e-6,'保存壁面与体场不同态。');
    % 模块四：继续研究使用独立、受进程锁定的谱系；局部步号绝不能冒充累计步号。
    continuation=fullfile(folder,'continuation.json');
    if isfile(continuation)
        path=fullfile(receipt.request.workDirectory,'continuation.json');h=fileSHA(continuation);
        token=regexp(receiptText,[regexptranslate('escape',jsonencode(path)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
        assert(numel(token)==1 && strcmp(h,token{1}{1}),'继续计算身份未被实际进程锁定或已改变。');
        lineage=jsondecode(fileread(continuation));before=lineage.completedUpdatesBeforeRun;
        validateattributes(before,{'numeric'},{'scalar','finite','integer','positive'});
        assert(before==lineage.sourceCumulativeSavedIteration+1 && lineage.samePhysicsAndNumerics && lineage.sameMesh);
        audit.continuation=lineage;audit.completedUpdatesBeforeRun=before;
        audit.cumulativeCheckpointIteration=before+audit.checkpointIteration;
        audit.cumulativeTerminalIteration=before+audit.terminalIteration;
        audit.successfulEngineSecondsThroughRun=lineage.priorSuccessfulEngineSeconds+receipt.elapsedSeconds;
        audit.evidenceSHA256.continuation_json=h;
    end
    % 改网格初态只记录来源，不能把新网格局部步号接成同网格累计步号。
    transfer=fullfile(folder,'initial_flow.transfer.json');
    if isfile(transfer)
        assert(~isfile(continuation),'改网格迁移与同网格续算身份冲突。');
        path=fullfile(receipt.request.workDirectory,'initial_flow.transfer.json');h=fileSHA(transfer);
        token=regexp(receiptText,[regexptranslate('escape',jsonencode(path)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
        assert(numel(token)==1 && strcmp(h,token{1}{1}),'迁移说明未被真实进程锁定。');
        lineage=jsondecode(fileread(transfer));
        assert(strcmp(lineage.kind,'normal-column-state-transfer') && ~lineage.sameMesh && lineage.initialGuessOnly);
        assert(strcmp(fileSHA(fullfile(folder,'initial_flow.csv')),lineage.outputSHA256),'迁移初态已改变。');
        audit.meshTransfer=lineage;audit.sourceCumulativeIteration=lineage.sourceCumulativeIteration;
        audit.evidenceSHA256.initial_flow_transfer_json=h;
    end
    % 新物性初猜不属于同物理续算；守恒能量零点改变后必须从局部步号记录。
    propertyTransfer=fullfile(folder,'coolprop_initialization.json');
    if isfile(propertyTransfer)
        assert(~isfile(continuation) && ~isfile(transfer),'物性转换与其他续算身份冲突。');
        path=fullfile(receipt.request.workDirectory,'coolprop_initialization.json');h=fileSHA(propertyTransfer);
        token=regexp(receiptText,[regexptranslate('escape',jsonencode(path)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
        assert(numel(token)==1 && strcmp(h,token{1}{1}),'物性转换未绑定真实进程。');
        lineage=jsondecode(fileread(propertyTransfer));
        assert(strcmp(lineage.kind,'same-mesh-primitive-to-HEOS-Air-initial-guess'));
        assert(lineage.sameMesh && ~lineage.samePhysicsAndNumerics && ~lineage.strictlyConservativeRemap);
        assert(lineage.newPhysicsLocalIterationStartsAt==0 && max(cell2mat(struct2cell(lineage.errors)))<1e-8);
        assert(strcmp(fileSHA(fullfile(folder,'initial_flow.csv')),lineage.initialStateSHA256));
        audit.propertyTransfer=lineage;audit.sourceCumulativeIteration=lineage.sourceCumulativeIteration;
        audit.evidenceSHA256.coolprop_initialization_json=h;
    end
    % SST转换只记录源模型的步号，绝不能把层流步号加到新模型历史。
    if sst && isfile(fullfile(folder,'sst_initialization.json'))
        path=fullfile(receipt.request.workDirectory,'sst_initialization.json');h=fileSHA(path);
        token=regexp(receiptText,[regexptranslate('escape',jsonencode(path)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
        assert(numel(token)==1 && strcmp(h,token{1}{1}),'SST初态来源未绑定实际进程。');
        assert(~isfile(continuation) && ~isfile(transfer) && ~isfile(propertyTransfer));
        lineage=jsondecode(fileread(path));assert(lineage.newPhysicsLocalIterationStartsAt==0 && ~lineage.samePhysicsAndNumerics);
        audit.turbulenceTransfer=lineage;audit.sourceModelLocalIteration=lineage.sourceLocalIteration;
    end
    audit.checkpointIdentityConfirmed=true;audit.status='checkpoint-and-history-audited';
    audit.reason='保存迭代已绑定且重启动EOS一致；残差通过与物理资格仍分别判断。';
catch ex
    audit.status='evidence-rejected';audit.reason=ex.message;
    audit.checkpointIdentityConfirmed=false;audit.checkpointResidualPassed=false;audit.terminalResidualPassed=false;
end
end

function value=number(text,key)
token=regexp(text,['(?m)^\s*' key '\s*=\s*([-+0-9.eE]+)'],'tokens','once');assert(~isempty(token),['缺少' key]);
value=str2double(token{1});assert(isscalar(value) && isfinite(value));
end

function h=fileSHA(file)
% SHA绑定文件身份，不改写原结果，也不扫描历史输出目录。
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),b=fread(fid,1024*1024,'*uint8');md.update(b);end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
