function receipt=exportWallPairEvidence(data,path)
%EXPORTWALLPAIREVIDENCE 从当前单时刻工作台账本导出同一原面上的壁压和热流。
% 模块一：只读取当前证据页已经建立的账本，确认两列场齐备且各自属于单一时刻。
ledger=data.evidenceLedger;
keys={ledger.fields.key};
p=ledger.fields(strcmp(keys,'pressure'));
q=ledger.fields(strcmp(keys,'heatFlux'));
assert(isscalar(p)&&isscalar(q),'RAH:WallPair:Fields','当前快照须同时含壁压和热流字段。');
assert(~(isfield(p,'history')&&~isempty(p.history)) && ~(isfield(q,'history')&&~isempty(q.history)), ...
    'RAH:WallPair:History','壁压与热流需要来自同一单时刻快照。');
n=ledger.faceCount;
assert(n==size(data.faces,1)&&n==numel(data.faceArea),'RAH:WallPair:Faces','原面数量不一致。');

% 模块二：原数据的单位和符号直接沿用。给定、超域、研究计算和工程计算分列。
t=table((1:n)',string(data.region(:)),compose("%.17g",data.faceArea(:)), ...
    'VariableNames',{'face_id','region','area_m2'});
t=addField(t,p,'pressure',ledger);
t=addField(t,q,'heat_flux',ledger);
t.joint_engineering=(p.availability(:)==1 & q.availability(:)==1);
t.joint_verified=(ismember(p.availability(:),[1 2]) & ismember(q.availability(:),[1 2]) & ...
    p.qualification(:)==2 & q.qualification(:)==2);
assert(isequal(find(t.joint_engineering)-1,ledger.joint.pressureHeatFluxEngineeringFaceIDs(:)) && ...
       isequal(find(t.joint_verified)-1,ledger.joint.pressureHeatFluxQualifiedFaceIDs(:)), ...
       'RAH:WallPair:Ledger','逐面联合状态与当前账本不一致。');
% 每行保留工况和源身份，筛选或截取 CSV 后仍能独立辨认来源。
caption='';if isfield(data.resultCondition,'caption'),caption=data.resultCondition.caption;end
t.condition=repmat(string(caption),n,1);
t.geometry_sha256=repmat(string(ledger.identity.geometrySHA256),n,1);
t.condition_sha256=repmat(string(ledger.identity.conditionSHA256),n,1);
t.source_sha256=repmat(string(ledger.identity.sourceSHA256),n,1);
t.source_traceable=repmat(ledger.identity.sourceTraceable,n,1);
t.synthetic=repmat(ledger.identity.synthetic,n,1);
% 热敏珠换算表给出原生参照及比例，供取出单行后复核压力与热流的来源。
if isfield(data,'thermistorLead') && isfield(data.thermistorLead,'thermal')
    thermal=data.thermistorLead.thermal;
    t.reference_frame_sha256=repmat(string(thermal.nativeFrameSHA256),n,1);
    t.reference_wall_K=repmat(thermal.referenceWallTemperatureK,n,1);
    t.temperature_difference_ratio=repmat(thermal.scale,n,1);
elseif isfield(data,'thermistorWire')
    % 非均温对照逐行记录源帧、边界输入、原面温度和两种边界的热流差。
    w=data.thermistorWire;tw=data.fields(strcmp({data.fields.key},'wallInput'));
    dq=data.fields(strcmp({data.fields.key},'wireHeatFluxChange'));
    t.source_frame_sha256=repmat(string(w.sourceFrameSHA256),n,1);
    t.uniform_reference_frame_sha256=repmat(string(w.baseFrame),n,1);
    t.wall_boundary_sha256=repmat(string(w.boundarySHA256),n,1);
    t.wall_input_K=tw.values(:);t.change_from_uniform_kW_m2=dq.values(:);
    if isfield(w,'thermometry')
        m=w.thermometry;t.resistance_temperature_K=repmat(m.internalK,n,1);t.inferred_bead_surface_K=repmat(m.surfaceK,n,1);
        t.uniform_reference_temperature_ratio=repmat(m.scale,n,1);t.source_field_temperature_ratio=repmat(m.fieldFactor,n,1);
        t.thermometry_reference_frame_sha256=repmat(string(m.referenceFrame),n,1);
        change=data.fields(strcmp({data.fields.key},'thermometryHeatFluxChange'));t.thermometry_change_kW_m2=change.values(:);
    end
    t.experiment_doi=repmat(string(w.experiment.doi),n,1);
elseif isfield(data,'measuredPoint')
    % 截取单行时仍可辨认当前 12 K 场量路径，以及论文实际测量的空间范围。
    m=data.measuredPoint;ref=m.scaledFrom;ratio=m.scaledFrom.temperatureRatio;
    if strcmp(m.action,'archive'),ref=m.nativeCheck;ratio=1;end
    if isfield(m,'coldSource'),ref=m.coldSource;ratio=1;end
    t.reference_frame_sha256=repmat(string(ref.frame),n,1);
    t.reference_wall_K=repmat(ref.wallK,n,1);t.temperature_difference_ratio=repmat(ratio,n,1);
    t.experiment_doi=repmat(string(m.experiment.doi),n,1);
    t.measurement_scope=repmat("图 7 珠部平均换热；逐面实验精度未知",n,1);
end

% 模块三：先写同目录临时文件，完整写出后再替换目标文件。
[folder,name,ext]=fileparts(char(path));if isempty(folder),folder=pwd;end
if isempty(ext),ext='.csv';end
assert(isfolder(folder)&&~isempty(name)&&strcmpi(ext,'.csv'),'RAH:WallPair:Path','请选择现有目录中的 CSV 文件。');
target=fullfile(folder,[name ext]);assert(~isfolder(target),'RAH:WallPair:Path','导出目标是目录。');
tmp=[tempname(folder) '.csv'];
clean=onCleanup(@()discardTemporary(tmp)); %#ok<NASGU>
writetable(t,tmp,'Encoding','UTF-8','QuoteStrings','all');
[ok,msg]=movefile(tmp,target,'f');assert(ok,'RAH:WallPair:Write','无法保存原面表。%s',msg);
receipt=struct('csvPath',target,'rows',n,'joint',ledger.joint,'identity',ledger.identity);
end

function t=addField(t,f,prefix,ledger)
% 使用原场量的算法字典与原始缺项说明。有限的超域值仍须保留其状态。
value=compose("%.17g",f.values(:));value(~isfinite(f.values(:)))="";
t.(prefix)=value;
t.([prefix '_unit'])=repmat(string(f.unit),height(t),1);
availability=string(ledger.availabilityNames(:));qualification=string(ledger.qualificationNames(:));
t.([prefix '_availability'])=availability(1+double(f.availability(:)));
t.([prefix '_qualification'])=qualification(1+double(f.qualification(:)));
if isfield(f,'methodNames')&&~isempty(f.methodNames)&&isfield(f,'methodIndex')&&~isempty(f.methodIndex)
    methods=string(f.methodNames(:));t.([prefix '_method'])=methods(1+double(f.methodIndex(:)));
else
    t.([prefix '_method'])=repmat("来源未记录",height(t),1);
end
reason=string(f.reason(:));
if isfield(f,'missingReason')&&~isempty(f.missingReason)
    original=string(f.missingReason(:));use=(f.availability(:)==0|f.availability(:)==4)&strlength(original)>0;
    reason(use)=original(use);
end
t.([prefix '_reason'])=reason;
end

function discardTemporary(path)
% 清理本次导出创建的临时文件。
if isfile(path),delete(path);end
end
