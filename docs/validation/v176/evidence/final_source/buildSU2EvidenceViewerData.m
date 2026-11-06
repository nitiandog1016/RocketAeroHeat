function data=buildSU2EvidenceViewerData(result,evidenceRecords)
%BUILDSU2EVIDENCEVIEWERDATA 把已保存SU2原面结果接入独立区域证据页。
% 模块一：只读内存快照，不读取体场、不调用映射/求解器、不修改原值。
% 模块二：原生p/q/流体温度/保守热负荷各保留身份；输入Tw不是材料温度解。
% 模块三：只有既存单位与帧身份审核齐备才标SI，旧未审核结果显示raw。
% 模块四：几何标签只能显式传入，有限研究值不自动授予物理资格。
if nargin<2
    if isstruct(result)&&isscalar(result)&&isfield(result,'evidenceRecords'),evidenceRecords=result.evidenceRecords;
    else,evidenceRecords=[];end
end
must(isstruct(result)&&isscalar(result),'需要一个已保存的SU2原面结果。');
must(isfield(result,'V')&&isfield(result,'F'),'结果必须保留原V/F。');
V=double(result.V);F=double(result.F);n=size(F,1);
must(size(V,2)==3&&size(F,2)==3&&all(isfinite(V(:)))&&all(isfinite(F(:)))&& ...
    all(F(:)==fix(F(:)))&&all(F(:)>=1&F(:)<=size(V,1)),'原V/F格式不合法。');
normal=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
area=sqrt(sum(normal.^2,2))/2;
if isfield(result,'faceArea'),area=double(result.faceArea(:));end
data=struct('kind','aerothermal-model','schemaVersion',1,'vertices',V,'faces',F-1, ...
    'faceArea',area,'fields',struct([]),'coordinateNote','原模型SI坐标；原节点、面元顺序与比例保持', ...
    'caseNote','既存原面研究快照；有限值不等于工程适用或物理验证通过');
data.region=repmat({'未分类'},n,1);
data.regionSource='原研究结果没有区域标签；不按坐标推造分离、再附、底面或突起区域';
if isfield(result,'geometryRegion')
    data.region=cellstr(string(result.geometryRegion(:)));
    must(numel(data.region)==n,'geometryRegion必须逐原面记录。');
    data.regionSource='显式原面几何标签；只供区域查看，不转移流动或热资格';
    if isfield(result,'geometryRegionSource'),data.regionSource=char(result.geometryRegionSource);end
end

[siVerified,sourceFrame]=auditedSource(result);
data.unitsVerified=siVerified;
if siVerified
    pressureUnit='Pa';heatUnit='W/m²';temperatureUnit='K';
    data.unitNote='沿用既存verifiedCaseReview、当前checkpoint与热状态接口审核的SI单位；本页未重读体场';
else
    pressureUnit='raw pressure';heatUnit='raw heat flux';temperatureUnit='raw temperature';
    data.unitNote='单位审核合同不完整；保留原始数值与raw单位，不猜测SI，不授物理资格';
end
solver=textOr(result,'solverMethod','源快照未记录求解方法');
projection=textOr(result,'projectionMethod','源快照未记录原面映射方法');
researchScope=['既存原面研究结果；' projection '；残差、近壁分辨率、物性及实验适用性须独立审核'];
add('pressure','表面压力',pressureUnit,values('pressureRaw'),1,solver,researchScope, ...
    '原面压力研究值；不以热流覆盖率或接口通过授予压力物理资格');
qInput=isTrue(result,'heatFluxIsBoundaryInput');
if qInput,qState=2;qMethod='热边界给定热流';qScope='给定热流边界输入；不是独立预测的对流热流';
else,qState=1;qMethod=solver;qScope=researchScope;end
add('heatFlux','对流热流密度',heatUnit,values('heatFluxRaw'),qState,qMethod,qScope, ...
    '原有符号热流保持；零值和负值不作为缺失值，物理资格独立审核');
add('fluidTemperature','壁邻流体节点平均温度',temperatureUnit,values('fluidNodeTemperatureMeanRaw'),1, ...
    [solver '；' projection],researchScope,textOr(result,'temperatureIdentity', ...
    '原面平均流体节点状态；不是给定Tw，也不是材料壁温解'));
add('edgeTemperature','外缘温度（未提取）',temperatureUnit,nan(n,1),0, ...
    '未提取适用的局部外缘状态','当前结果未提供经过适用性审核的逐面外缘温度', ...
    '分离、再附或回流区可能没有唯一局部外缘定义；壁邻流体温度不能代填外缘温度');
add('recoveryTemperature','恢复温度（未提取）',temperatureUnit,nan(n,1),0, ...
    '未提取适用的局部恢复状态','当前结果未提供经过适用性审核的逐面恢复温度', ...
    '需要适用的局部外缘、恢复关系及热物性；回流区可能无唯一定义，不由壁邻流体温度补齐');
add('heatLoadDensity','保守热负荷面密度',heatUnit,values('conservativeHeatLoadDensityRaw'),1, ...
    textOr(result,'heatLoadMethod','源快照未记录保守热负荷方法'),researchScope, ...
    '原生节点合成面积对应的保守热负荷；与原面平均热流是不同量，不能互相替代');
add('wallInput','热边界输入壁温',temperatureUnit,values('wallTemperatureInput'),2, ...
    '源结果明确保存的热边界输入', '给定或耦合传入的热边界输入；不是本次材料壁温求解', ...
    '输入面积100%不表示求解了壁体温度场；本页不由流体温度或初始温度补值');
add('wallSolved','求解壁温',temperatureUnit,nan(n,1),0,'未提供独立材料壁温解', ...
    '当前流体原面结果未提供独立材料热响应状态及来源', ...
    '流体节点温度和热边界输入均不能替代材料壁温解');

condition=conditionSnapshot(result,values('wallTemperatureInput'));
physical=condition;physical=rmfield(physical,{'caption','source','wallIdentity'});
if isfield(result,'thermalBoundary')&&isstruct(result.thermalBoundary)
    boundary=result.thermalBoundary;
    for key={'mode','specification','heatFluxIsBoundaryInput'}
        if isfield(boundary,key{1}),physical.thermalBoundary.(key{1})=boundary.(key{1});end
    end
end
if isfield(result,'receipt')&&isstruct(result.receipt)
    receipt=result.receipt;
    if isfield(receipt,'physicsModel'),physical.physicsModel=receipt.physicsModel;end
    if isfield(receipt,'thermodynamicModel')&&isstruct(receipt.thermodynamicModel)
        model=receipt.thermodynamicModel;
        for key={'thermodynamicModel','gasConstant','temperatureRangeK','energyEntropyReference','turbulent','flowFlux','conductivity'}
            if isfield(model,key{1}),physical.thermodynamicModel.(key{1})=model.(key{1});end
        end
    end
end
data.resultCondition=condition;data.caseNote=[condition.caption '。' data.unitNote];
data.evidenceIdentity=struct('origin','research','sourceTraceable',siVerified, ...
    'synthetic',false,'sourceFrame',sourceFrame,'physicalCondition',physical);
data.evidenceRecords=evidenceRecords;
data.evidenceLedger=buildAerothermalEvidenceLedger(data,evidenceRecords);
data.totalArea=data.evidenceLedger.totalArea;

    function value=values(name)
        value=nan(n,1);
        if isfield(result,name),value=double(result.(name)(:));end
        % 旧研究入口允许明确给定的均匀壁温标量；仅此边界输入可广播。
        if strcmp(name,'wallTemperatureInput')&&isscalar(value),value=repmat(value,n,1);end
        must(numel(value)==n&&isreal(value)&&~any(isinf(value)),['原字段 ' name ' 必须逐面保持有限值或NaN。']);
    end
    function add(key,label,unit,value,finiteState,method,scope,note)
        state=finiteState*double(isfinite(value));
        reasons=repmat({''},n,1);reasons(state==0)={'源快照未提供此变量的有限原面结果；不以其他变量补齐'};
        if strcmp(key,'wallSolved'),reasons(:)={'当前流体结果没有独立材料壁温解；未执行壁体温度求解'};end
        outside=false(n,1);
        if isfield(result,'evidenceOutOfDomainFace')&&isstruct(result.evidenceOutOfDomainFace)&& ...
                isfield(result.evidenceOutOfDomainFace,key)
            outside=result.evidenceOutOfDomainFace.(key);
        end
        entry=struct('key',key,'label',label,'unit',unit,'values',value,'status',state, ...
            'methodNames',{{method}},'methodIndex',zeros(n,1),'scope',scope, ...
            'scopeNames',{{scope}},'scopeIndex',zeros(n,1),'missingReason',{reasons}, ...
            'stateNote',note,'evidenceOrigin','research','outOfDomainFace',outside);
        if isempty(data.fields),data.fields=entry;else,data.fields(end+1)=entry;end
    end
end

function [verified,frame]=auditedSource(result)
% 此处绑定已保存审核身份，不把本次只读展示称为原始文件重验。
frame=struct('kind','saved-su2-original-surface','auditReperformedByAdapter',false);
for key={'resultDirectory','sourceRevision','solveId','solverMethod','projectionMethod', ...
        'actualRunControls','sourceCheckpointIdentity'}
    if isfield(result,key{1}),frame.(key{1})=result.(key{1});end
end
verified=false;
if ~isfield(result,'verifiedCaseReview')||~isstruct(result.verifiedCaseReview)|| ...
        ~isfield(result,'convergenceAudit')||~isstruct(result.convergenceAudit)|| ...
        ~isfield(result,'thermalStateAudit')||~isstruct(result.thermalStateAudit),return;end
review=result.verifiedCaseReview;audit=result.convergenceAudit;thermal=result.thermalStateAudit;
frame.verifiedCaseReview=review;
if isfield(audit,'evidenceSHA256'),frame.currentCheckpointEvidenceSHA256=audit.evidenceSHA256;end
for key={'checkpointIteration','completedUpdatesBeforeRun','cumulativeCheckpointIteration','nativeCheckpointIdentity'}
    if isfield(audit,key{1}),frame.(key{1})=audit.(key{1});end
end
for key={'auditSHA256','wallCSV_SHA256','checkpointLocalIteration','flowPrimitivePhase','turbulencePhase'}
    if isfield(thermal,key{1}),frame.thermalState.(key{1})=thermal.(key{1});end
end
if ~isTrue(audit,'checkpointIdentityConfirmed')||~isTrue(thermal,'interfacePassed')|| ...
        ~hashField(review,'manifestSHA256')||~hashField(review,'meshSHA256')|| ...
        ~hashField(thermal,'auditSHA256')||~isfield(audit,'evidenceSHA256')||~isstruct(audit.evidenceSHA256),return;end
evidence=audit.evidenceSHA256;
for key={'engine_process_json','case_cfg','surface_flow_vtk','restart_flow_csv'}
    if ~hashField(evidence,key{1}),return;end
end
verified=true;
end

function condition=conditionSnapshot(result,wall)
condition=struct('source','结果自己的冻结工况；不读取待编辑配置','Mach',NaN,'alphaDeg',NaN, ...
    'betaDeg',NaN,'altitudeM',NaN,'freestreamPressurePa',NaN,'freestreamTemperatureK',NaN, ...
    'wallInputType','unknown','wallUniformValueK',NaN,'wallMinK',NaN,'wallMaxK',NaN, ...
    'wallKnownFaces',nnz(isfinite(wall)),'wallUnknownFaces',nnz(~isfinite(wall)), ...
    'wallIdentity','热边界输入；不是材料壁温解');
if isfield(result,'flight')&&isstruct(result.flight)
    from={'Mach','alpha','beta','altitude','pressure','temperature'};
    to={'Mach','alphaDeg','betaDeg','altitudeM','freestreamPressurePa','freestreamTemperatureK'};
    for k=1:numel(from)
        if isfield(result.flight,from{k})
            value=result.flight.(from{k});
            if isnumeric(value)&&isscalar(value)&&isfinite(value),condition.(to{k})=double(value);end
        end
    end
end
known=isfinite(wall);
if any(known)
    condition.wallMinK=min(wall(known));condition.wallMaxK=max(wall(known));condition.wallInputType='perFace';
    if all(known)&&all(wall==wall(1)),condition.wallInputType='uniform';condition.wallUniformValueK=wall(1);end
end
condition.caption=sprintf('既存研究快照：M %s，α %s°，β %s°；p∞ %s Pa，T∞ %s K；输入壁温 %s，未知 %d 面', ...
    number(condition.Mach),number(condition.alphaDeg),number(condition.betaDeg), ...
    number(condition.freestreamPressurePa),number(condition.freestreamTemperatureK),condition.wallInputType,condition.wallUnknownFaces);
end

function value=textOr(s,key,fallback)
value=fallback;
if isfield(s,key)&&(ischar(s.(key))||(isstring(s.(key))&&isscalar(s.(key))))&&~isempty(s.(key)),value=char(s.(key));end
end
function yes=isTrue(s,key)
yes=isfield(s,key)&&(islogical(s.(key))||isnumeric(s.(key)))&&isscalar(s.(key))&&isequal(double(s.(key)),1);
end
function yes=hashField(s,key)
yes=isfield(s,key)&&(ischar(s.(key))||(isstring(s.(key))&&isscalar(s.(key))))&& ...
    ~isempty(regexp(char(s.(key)),'^[0-9a-fA-F]{64}$','once'));
end
function text=number(value)
if isfinite(value),text=sprintf('%.8g',value);else,text='未记录';end
end
function must(condition,message)
if ~condition,error('RocketAeroHeat:SU2Evidence:Contract','%s',message);end
end
