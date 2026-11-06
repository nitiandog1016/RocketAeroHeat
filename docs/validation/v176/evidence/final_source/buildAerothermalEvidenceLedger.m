function ledger=buildAerothermalEvidenceLedger(data,evidenceRecords)
%BUILDAEROTHERMALEVIDENCELEDGER 从既有逐面元快照构造可用性/物理资格两层台账。
% 模块一：核对原V/F、真实面积、字段与冻结物理条件；不调用求解器或旧viewer。
% 模块二：旧status与值原样保留；工程/研究/给定/未知/超域是独立可用性。
% 模块三：只有准确SHA小证据中的逐变量规则可授资格；失败优先，不按passed猜结论。
% 模块四：按原faceID并集汇总区域及p/q交集，不缩分母、不重复计数。
if nargin<2
    if isstruct(data)&&isscalar(data)&&isfield(data,'evidenceRecords'),evidenceRecords=data.evidenceRecords;
    else,evidenceRecords=[];end
end
must(isstruct(data)&&isscalar(data),'需要一个canonical viewerData。');
for key={'vertices','faces','faceArea','region','fields','resultCondition'}
    must(isfield(data,key{1}),['缺少数据字段 ' key{1}]);
end
V=double(data.vertices);F=double(data.faces);area=double(data.faceArea(:));n=size(F,1);
must(ismatrix(V)&&size(V,2)==3&&~isempty(V)&&all(isfinite(V(:))),'原节点必须为有限Nx3。');
must(ismatrix(F)&&size(F,2)==3&&n>0&&all(isfinite(F(:)))&& ...
    all(F(:)==fix(F(:)))&&all(F(:)>=0&F(:)<size(V,1)),'原三角面必须使用零起点合法节点编号。');
crosses=cross(V(F(:,2)+1,:)-V(F(:,1)+1,:),V(F(:,3)+1,:)-V(F(:,1)+1,:),2);
exactArea=sqrt(sum(crosses.^2,2))/2;
must(numel(area)==n&&all(isfinite(area)&area>0)&&all(exactArea>0)&& ...
    max(abs(area-exactArea)./max(exactArea,realmin))<2e-10,'faceArea与全部原三角面积不一致。');
region=cellstr(string(data.region(:)));must(numel(region)==n,'区域标签数量与原面不符。');
must(isstruct(data.fields)&&~isempty(data.fields),'至少需要一个现有字段。');
fields=data.fields;keys={fields.key};must(numel(unique(keys))==numel(keys),'变量键不能重复。');
origin='engineering';sourceFrame=struct('kind','unrecorded-memory-snapshot');traceable=false;synthetic=false;
if isfield(data,'evidenceIdentity')
    identityInput=data.evidenceIdentity;must(isstruct(identityInput)&&isscalar(identityInput),'evidenceIdentity必须为结构。');
    if isfield(identityInput,'origin'),origin=char(identityInput.origin);end
    if isfield(identityInput,'sourceFrame'),sourceFrame=identityInput.sourceFrame;end
    if isfield(identityInput,'sourceTraceable'),traceable=isequal(identityInput.sourceTraceable,true);end
    if isfield(identityInput,'synthetic'),synthetic=isequal(identityInput.synthetic,true);end
else,identityInput=struct;end
must(ismember(origin,{'engineering','research'}),'未知结果来源类别。');
physical=stripDisplay(data.resultCondition);
if isfield(identityInput,'physicalCondition'),physical=stripDisplay(identityInput.physicalCondition);end
boundaryInputs=struct;
for k=1:numel(fields)
    f=fields(k);must(ischar(f.key)||isstring(f.key),'变量键必须为文本。');
    must(isvarname(char(f.key)),'变量键必须可作为MATLAB结构字段。');
    for name={'label','unit','values','status'}
        must(isfield(f,name{1}),['字段 ' char(f.key) ' 缺少 ' name{1}]);
    end
    values=double(f.values(:));status=double(f.status(:));
    must(numel(values)==n&&isreal(values)&&~any(isinf(values)),'字段必须逐原面保存有限值或NaN。');
    must(numel(status)==n&&all(ismember(status,[0 1 2])),'旧字段状态只接受0/1/2。');
    must(all(isnan(values(status==0)))&&all(isfinite(values(status~=0))),'旧状态与有限/NaN值不一致。');
    if any(status==2)
        boundaryInputs.(char(f.key))=struct('faceIDs',find(status==2)-1,'values',values(status==2),'unit',char(f.unit));
    end
end
identity=struct('geometrySHA256',hashValue(struct('vertices',V,'faces',F)), ...
    'conditionSHA256',hashValue(struct('physical',physical,'prescribedInputs',boundaryInputs)), ...
    'sourceSHA256',hashValue(struct('origin',origin,'sourceFrame',sourceFrame)), ...
    'sourceTraceable',traceable,'origin',origin,'synthetic',synthetic);
for k=1:numel(fields)
    f=fields(k);values=double(f.values(:));status=double(f.status(:));
    fieldOrigin=origin;
    if isfield(f,'evidenceOrigin')&&~isempty(f.evidenceOrigin),fieldOrigin=char(f.evidenceOrigin);end
    must(ismember(fieldOrigin,{'engineering','research'}),'未知逐变量来源类别。');
    availability=zeros(n,1);availability(status==2)=3;
    if strcmp(fieldOrigin,'research'),availability(status==1)=2;else,availability(status==1)=1;end
    outside=false(n,1);
    if isfield(f,'outOfDomainFace')&&~isempty(f.outOfDomainFace)
        v=f.outOfDomainFace(:);must(numel(v)==n&&all(ismember(v,[0 1])),'超域掩码不合法。');outside=logical(v);
        % 明确超域也可能没有可用数值；保留NaN和旧status=0，不能为了分类补值。
    end
    availability(outside)=4;qualification=zeros(n,1);qualification(outside)=1;
    reason=repmat({'缺少与本变量、原面、工况和来源匹配的独立物理证据'},n,1);
    reason(availability==0)={'本变量未知；不能由其他变量或零值补齐'};
    reason(availability==3)={'给定边界输入，不是本变量的独立求解或物理验证'};
    reason(outside)={'原结果明确超出适用域；保留实际值或NaN，物理资格失败'};
    fields(k).legacyStatus=f.status;fields(k).availability=availability;
    fields(k).qualification=qualification;fields(k).reason=reason;
    fields(k).evidenceIDs=repmat({cell(0,1)},n,1);
    % valueSHA同时绑定数值身份与实际算法来源；等值换方法也必须重新取证。
    fieldBinding=struct('key',char(f.key),'unit',char(f.unit),'values',values, ...
        'legacyStatus',status,'evidenceOrigin',fieldOrigin,'outOfDomainFace',outside);
    for name={'methodNames','methodIndex','scope','scopeNames','scopeIndex'}
        if isfield(f,name{1}),fieldBinding.(name{1})=f.(name{1});end
    end
    fields(k).valueSHA256=hashValue(fieldBinding);
end

% 先对记录身份去重；同id冲突时整组拒绝，避免依赖输入顺序留下先到的通过。
if isempty(evidenceRecords),records={};
elseif isstruct(evidenceRecords),records=num2cell(evidenceRecords);
elseif iscell(evidenceRecords),records=evidenceRecords;
else,error('RocketAeroHeat:Evidence:Contract','evidenceRecords必须为结构数组或cell。');end
records=records(:);ids=cell(size(records));fingerprints=cell(size(records));preErrors=cell(size(records));
for i=1:numel(records)
    try
        r=records{i};must(isstruct(r)&&isscalar(r),'记录必须为单结构。');
        must(strcmp(char(r.schema),'RAH.AerothermalEvidenceRecord.v1'),'证据记录schema不支持。');
        ids{i}=textValue(r.id,'id');must(~isempty(ids{i}),'证据id不能为空。');
        field=textValue(r.variable,'variable');faces=faceSet(r.faceIDs,n);
        must(isfield(r,'bindings')&&isstruct(r.bindings),'记录缺少bindings。');
        fingerprint=struct('id',ids{i},'variable',field,'faceIDs',faces, ...
            'bindings',r.bindings,'proofSHA256',lower(textValue(r.proofSHA256,'proofSHA256')));
        fingerprints{i}=hashValue(fingerprint);
    catch ex,preErrors{i}=ex.message;end
end
accepted=struct([]);rejected=struct([]);duplicates=struct([]);
for i=1:numel(records)
    r=records{i};id=ids{i};if isempty(id),id=sprintf('invalid-record-%d',i);end
    try
        must(isempty(preErrors{i}),preErrors{i});
        group=find(strcmp(ids,ids{i}));
        must(numel(unique(fingerprints(group)))==1,'相同证据id包含冲突内容；整组拒绝。');
        if i~=group(1)
            duplicates=append(duplicates,struct('id',id,'recordIndex',i));continue;
        end
        variable=char(r.variable);k=find(strcmp(keys,variable));must(isscalar(k),'证据变量未在当前快照中出现。');
        faces=faceSet(r.faceIDs,n);ix=faces+1;f=fields(k);
        binding=struct('geometrySHA256',identity.geometrySHA256,'conditionSHA256',identity.conditionSHA256, ...
            'sourceSHA256',identity.sourceSHA256,'valueSHA256',f.valueSHA256);
        must(sameBindings(r.bindings,binding),'证据几何/工况/源帧/变量值SHA不匹配。');
        must(traceable,'当前源身份未被明确绑定；内存快照不能自动取得物理资格。');
        must(all(ismember(f.availability(ix),[1 2 4])),'证据不能把给定或未知值升级为计算物理资格。');
        path=textValue(r.proofFile,'proofFile');file=java.io.File(path);
        must(file.isAbsolute()&&isfile(path),'proofFile必须是实际存在的绝对路径小JSON。');
        sizeInfo=dir(path);must(isscalar(sizeInfo)&&sizeInfo.bytes>0&&sizeInfo.bytes<=1024^2,'证据JSON超过1MiB或为空。');
        proofSHA=fileHash(path);must(validHash(r.proofSHA256)&&strcmpi(proofSHA,char(r.proofSHA256)),'证据文件SHA不匹配。');
        proof=jsondecode(fileread(path));
        must(isstruct(proof)&&isscalar(proof)&&strcmp(char(proof.schema),'RAH.AerothermalQualificationProof.v1'),'物理证据schema不支持。');
        must(strcmp(textValue(proof.evidenceID,'evidenceID'),id)&&strcmp(char(proof.variable),variable), ...
            '证据文件的id或变量不匹配。');
        must(isequal(faceSet(proof.faceIDs,n),faces)&&sameBindings(proof.bindings,binding),'证据文件原面集或绑定不匹配。');
        basis=textValue(proof.basis,'basis');
        must(ismember(basis,{'matched-experiment','validated-reference','synthetic-fixture'}),'证据不是声明的物理比较类型。');
        must(~strcmp(basis,'synthetic-fixture')||synthetic,'合成证据只允许显式合成数据，不授真实场资格。');
        applicability=textValue(proof.applicability,'applicability');must(~isempty(strtrim(applicability)),'缺少适用范围。');
        must(isstruct(proof.reference)&&isscalar(proof.reference),'缺少独立参照来源。');
        reference=proof.reference;referenceType=textValue(reference.type,'reference.type');
        must(~isempty(referenceType)&&~isempty(textValue(reference.uri,'reference.uri'))&&validHash(reference.sha256), ...
            '独立参照需要明确类型、URI和SHA。');
        [status,criteria]=criterionDecision(proof.criteria,variable);
        label={'未审核','失败','通过'};
        for j=ix'
            fields(k).evidenceIDs{j}{end+1,1}=id;
            previous=fields(k).qualification(j);
            if status==1,fields(k).qualification(j)=1;
            elseif status==2&&previous~=1,fields(k).qualification(j)=2;end
            note=sprintf('证据%s：%s；%s',id,label{status+1},applicability);
            if isempty(fields(k).evidenceIDs{j})||numel(fields(k).evidenceIDs{j})==1
                if previous==1,fields(k).reason{j}=[fields(k).reason{j} '；' note];else,fields(k).reason{j}=note;end
            else,fields(k).reason{j}=[fields(k).reason{j} '；' note];end
        end
        accepted=append(accepted,struct('id',id,'variable',variable,'faceIDs',faces,'qualification',status, ...
            'basis',basis,'applicability',applicability,'proofFile',path,'proofSHA256',proofSHA, ...
            'reference',reference,'criteria',criteria,'bindings',binding));
    catch ex
        rejected=append(rejected,struct('id',id,'recordIndex',i,'reason',ex.message));
    end
end

% 面积用原faceID集合；同面多份证据不重复计面积，失败覆盖通过。
for k=1:numel(fields),fields(k).areaSummary=areas(area,fields(k).availability,fields(k).qualification);end
names=unique(region,'stable');regions=struct([]);
for j=1:numel(names)
    selected=find(strcmp(region,names{j}));summary=struct;
    for k=1:numel(fields)
        summary.(char(fields(k).key))=areas(area(selected),fields(k).availability(selected),fields(k).qualification(selected));
    end
    regions=append(regions,struct('name',names{j},'faceIDs',selected-1,'area',sum(area(selected)), ...
        'fields',summary,'joint',jointAreas(fields,area,selected)));
end
ledger=struct;ledger.schemaVersion=1;ledger.kind='aerothermal-evidence-ledger';ledger.identity=identity;
ledger.vertices=data.vertices;ledger.faces=data.faces;ledger.faceArea=data.faceArea;ledger.region=data.region;
ledger.resultCondition=data.resultCondition;ledger.sourceFrame=sourceFrame;ledger.totalArea=sum(area);ledger.faceCount=n;
ledger.fields=fields;ledger.regions=regions;ledger.joint=jointAreas(fields,area,(1:n)');
ledger.availabilityNames={'unknown','engineering','research','prescribed','out-of-domain'};
ledger.qualificationNames={'unknown','failed','passed'};
ledger.acceptedEvidence=accepted;ledger.rejectedEvidence=rejected;ledger.duplicateEvidence=duplicates;
ledger.acceptedEvidenceCount=numel(accepted);ledger.rejectedEvidenceCount=numel(rejected);ledger.duplicateEvidenceCount=numel(duplicates);
ledger.qualificationIsIndependentOfAvailability=true;ledger.legacyValuesAndStatusPreserved=true;
ledger.physicalQualificationIsFromExplicitBoundProofOnly=true;ledger.synthetic=synthetic;
ledger.scope='原面逐变量两层证据台账；统计不求解、不补值、不自授资格；合成数据通过不属于真实物理认证';
end

function value=areas(area,a,q)
total=sum(area);value=struct('totalArea',total,'engineeringArea',sum(area(a==1)), ...
    'researchArea',sum(area(a==2)),'prescribedArea',sum(area(a==3)),'unknownArea',sum(area(a==0)), ...
    'outOfDomainArea',sum(area(a==4)),'qualifiedArea',sum(area(q==2)), ...
    'qualificationFailedArea',sum(area(q==1)),'qualificationUnknownArea',sum(area(q==0)));
names=fieldnames(value);
for i=2:numel(names),value.([names{i} 'Fraction'])=value.(names{i})/total;end
must(abs(value.engineeringArea+value.researchArea+value.prescribedArea+value.unknownArea+value.outOfDomainArea-total)<=2e-12*total,'可用面积未闭合。');
must(abs(value.qualifiedArea+value.qualificationFailedArea+value.qualificationUnknownArea-total)<=2e-12*total,'资格面积未闭合。');
end

function value=jointAreas(fields,area,selected)
keys={fields.key};ip=find(strcmp(keys,'pressure'));iq=find(strcmp(keys,'heatFlux'));n=numel(area);
available=false(n,1);qualified=available;engineering=available;
if isscalar(ip)&&isscalar(iq)
    p=fields(ip);q=fields(iq);
    available=ismember(p.availability,[1 2])&ismember(q.availability,[1 2]);
    engineering=p.availability==1&q.availability==1;
    qualified=p.qualification==2&q.qualification==2&available;
end
total=sum(area(selected));value=struct('totalArea',total);
for pair={'Available','Qualified','Engineering'}
    name=pair{1};if strcmp(name,'Available'),mask=available;elseif strcmp(name,'Qualified'),mask=qualified;else,mask=engineering;end
    ids=selected(mask(selected));prefix=['pressureHeatFlux' name];
    value.([prefix 'Area'])=sum(area(ids));value.([prefix 'Fraction'])=sum(area(ids))/total;value.([prefix 'FaceIDs'])=ids-1;
end
end

function [status,result]=criterionDecision(criteria,variable)
must(isstruct(criteria)&&~isempty(criteria),'证据必须包含数值验收规则，不接受泛型passed。');
required={'condition-match','variable-reference-comparison','numerical-accuracy','model-region-applicability'};
if ismember(variable,{'heatFlux','heatLoadDensity'}),required=[required,{'thermal-boundary-properties','near-wall-resolution'}];end
names={criteria.name};must(numel(unique(names))==numel(names)&&all(ismember(required,names)),'缺少本变量必需规则或规则重名。');
outcomes=zeros(numel(criteria),1);result=struct([]);
for i=1:numel(criteria)
    c=criteria(i);metric=textValue(c.metric,'metric');operator=textValue(c.operator,'operator');unit=textValue(c.unit,'unit');
    must(~isempty(metric)&&~isempty(unit),'每项规则须给出指标与单位。');
    threshold=c.threshold;must(isnumeric(threshold)&&isscalar(threshold)&&isreal(threshold)&&isfinite(threshold),'规则阈值必须有限。');
    observed=c.observed;must(isempty(observed)||(isnumeric(observed)&&isscalar(observed)&&isreal(observed)&&~isinf(observed)),'实测指标必须为标量或未知。');
    must(ismember(operator,{'<=','>=','<','>','==','abs<='}),'规则比较算子不支持。');
    if isempty(observed)||isnan(observed),outcomes(i)=0;
    else
        switch operator
            case '<=',ok=observed<=threshold;case '>=',ok=observed>=threshold;
            case '<',ok=observed<threshold;case '>',ok=observed>threshold;
            case '==',ok=observed==threshold;case 'abs<=',ok=abs(observed)<=threshold;
        end
        if ok,outcomes(i)=2;else,outcomes(i)=1;end
    end
    result=append(result,struct('name',char(c.name),'metric',metric,'observed',observed, ...
        'operator',operator,'threshold',threshold,'unit',unit,'qualification',outcomes(i)));
end
if any(outcomes==1),status=1;elseif all(outcomes==2),status=2;else,status=0;end
end

function faces=faceSet(value,n)
must(isnumeric(value)&&isreal(value)&&isvector(value)&&~isempty(value),'证据原面集合不能为空。');
faces=double(value(:));must(all(isfinite(faces))&&all(faces==fix(faces))&&all(faces>=0&faces<n)&&all(diff(faces)>0), ...
    '证据faceIDs必须为零起点、严格升序、无重复的原面编号。');
end

function yes=sameBindings(a,b)
names={'geometrySHA256','conditionSHA256','sourceSHA256','valueSHA256'};yes=isstruct(a)&&isscalar(a);
for i=1:numel(names)
    name=names{i};yes=yes&&isfield(a,name)&&validHash(a.(name))&&strcmpi(char(a.(name)),b.(name));
    if ~yes,return;end
end
end

function yes=validHash(value)
yes=(ischar(value)||(isstring(value)&&isscalar(value)))&&~isempty(regexp(char(value),'^[0-9a-fA-F]{64}$','once'));
end

function text=textValue(value,label)
must(ischar(value)||(isstring(value)&&isscalar(value)),[label '必须为文本']);text=char(value);
end

function value=stripDisplay(value)
if ~isstruct(value),return;end
display={'caption','source','label','note','description','wallIdentity','stateNote','scope'};
remove=intersect(fieldnames(value),display);if ~isempty(remove),value=rmfield(value,remove);end
for name=fieldnames(value)'
    key=name{1};for i=1:numel(value),value(i).(key)=stripDisplay(value(i).(key));end
end
end

function value=ordered(value)
if isstruct(value)
    value=orderfields(value);
    for name=fieldnames(value)'
        key=name{1};for i=1:numel(value),value(i).(key)=ordered(value(i).(key));end
    end
elseif iscell(value),for i=1:numel(value),value{i}=ordered(value{i});end
end
end

function h=hashValue(value)
% 结构字段排序后绑定精确double表示；显示字符串已从物理condition中单独剔除。
md=java.security.MessageDigest.getInstance('SHA-256');bytes=unicode2native(jsonencode(ordered(value)),'UTF-8');md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function h=fileHash(path)
fid=fopen(path,'rb');must(fid>=0,'无法读取小证据文件。');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,1024^2+1,'*uint8');must(numel(bytes)<=1024^2,'读取期间证据超过1MiB。');
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end

function must(condition,message)
if ~condition,error('RocketAeroHeat:Evidence:Contract','%s',message);end
end
