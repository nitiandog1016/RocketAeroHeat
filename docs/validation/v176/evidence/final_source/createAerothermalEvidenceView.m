function view=createAerothermalEvidenceView(host,data,varargin)
%CREATEAEROTHERMALEVIDENCEVIEW 在同一工作台查看原面元区域结果与验证证据。
% 模块一：仅消费已保存的逐面元结果及证据账本，不调用流场或热求解器。
% data 来自 buildAerothermalViewerData 或 buildSU2EvidenceViewerData。
% 可选第三参数为显式证据记录；未知、给定、研究值和验证资格分别显示。
assert(isvalid(host),'RocketAeroHeat:EvidenceView:Host','区域页宿主无效。');
if isa(host,'matlab.ui.container.GridLayout'),parent=host;
else,parent=uigridlayout(host,[1 1]);parent.Padding=[0 0 0 0];parent.RowSpacing=0;parent.ColumnSpacing=0;end
html=uihtml(parent,'HTMLSource',fullfile(fileparts(mfilename('fullpath')),'web','evidence.html'));
html.UserData=struct('payload',[],'revision',0,'ready',false,'rendered',[], ...
    'error',[],'viewState',struct,'sectionInfo',[]);
html.DataChangedFcn=@receive;
view=struct('html',html,'host',host,'update',@update,'getState',@getState,'getData',@getData, ...
    'selectField',@selectField,'selectRegion',@selectRegion,'inspectFace',@inspectFace, ...
    'setRange',@setRange,'setSectionX',@setSectionX,'getSection',@getSection, ...
    'exportTo',@exportTo,'exportSectionTo',@exportSectionTo);
update(data,varargin{:});

    function update(next,varargin)
        % 模块二：来源身份变化立即清除旧选择和回执；同一身份可保留固定色标。
        records=[];
        if ~isempty(varargin),records=varargin{1};elseif isfield(next,'evidenceRecords'),records=next.evidenceRecords;end
        % 外部附带的旧账本不能替代本次几何/工况/来源绑定；原证据逐项重新审核。
        next.evidenceRecords=records;next.evidenceLedger=buildAerothermalEvidenceLedger(next,records);
        assert(isfield(next,'fields')&&~isempty(next.fields)&&isfield(next,'faces'), ...
            'RocketAeroHeat:EvidenceView:Data','区域页需要原面元和逐变量结果。');
        s=html.UserData;same=~isempty(s.payload)&&isequaln(s.payload.data.evidenceLedger.identity,next.evidenceLedger.identity);
        if same,v=s.viewState;else,v=struct('fieldKey',next.fields(1).key,'regionName','', ...
                'selectedFace',-1,'position',[],'ranges',struct,'sectionX',mean([min(next.vertices(:,1)),max(next.vertices(:,1))]));end
        if ~any(strcmp({next.fields.key},v.fieldKey)),v.fieldKey=next.fields(1).key;end
        if ~isfield(v.ranges,v.fieldKey),v.ranges.(v.fieldKey)=initialRange(next.evidenceLedger.fields(strcmp({next.evidenceLedger.fields.key},v.fieldKey)));end
        % 每次数据更新清空原面元选择，避免旧索引被误用到新的结果快照。
        v.regionName='';v.selectedFace=-1;v.position=[];s.viewState=v;s.sectionInfo=[];
        s.revision=s.revision+1;s.rendered=[];s.error=[];
        s.payload=struct('kind','aerothermal-evidence-view','schemaVersion',1,'revision',s.revision,'data',next,'viewState',v);
        html.UserData=s;html.Data=s.payload;
    end
    function receive(src,~)
        message=src.Data;s=src.UserData;if ~isstruct(message)||~isfield(message,'kind'),return;end
        if strcmp(message.kind,'ready'),s.ready=true;src.UserData=s;src.Data=s.payload;return;end
        if ~isfield(message,'revision')||message.revision~=s.revision,return;end
        if any(strcmp(message.kind,{'evidence-state','evidence-rendered','evidence-export','evidence-section-export'}))
            for key={'fieldKey','regionName','selectedFace','position','ranges','sectionX'}
                if isfield(message,key{1}),s.viewState.(key{1})=message.(key{1});end
            end
            if isfield(message,'sectionInfo'),s.sectionInfo=message.sectionInfo;end
            s.payload.viewState=s.viewState;
            if strcmp(message.kind,'evidence-rendered'),s.rendered=message;s.error=[];end
            src.UserData=s;
        elseif strcmp(message.kind,'evidence-error'),s.error=message;src.UserData=s;return;
        else,return;end
        if strcmp(message.kind,'evidence-export')
            [name,folder]=uiputfile('*.csv','导出当前变量证据CSV','aerothermal_evidence.csv');
            if ~isequal(name,0),exportTo(fullfile(folder,name));end
        elseif strcmp(message.kind,'evidence-section-export')
            [name,folder]=uiputfile('*.csv','导出当前截面证据CSV','aerothermal_evidence_section.csv');
            if ~isequal(name,0),exportSectionTo(fullfile(folder,name));end
        end
    end
    function publish(v)
        s=html.UserData;s.revision=s.revision+1;s.viewState=v;s.sectionInfo=[];s.rendered=[];s.error=[];
        s.payload.revision=s.revision;s.payload.viewState=v;html.UserData=s;html.Data=s.payload;
    end
    function selectField(key)
        s=html.UserData;key=char(key);idx=find(strcmp({s.payload.data.evidenceLedger.fields.key},key),1);
        assert(~isempty(idx),'RocketAeroHeat:EvidenceView:Field','未知场量。');v=s.viewState;v.fieldKey=key;
        if ~isfield(v.ranges,key),v.ranges.(key)=initialRange(s.payload.data.evidenceLedger.fields(idx));end
        publish(v);
    end
    function selectRegion(name)
        s=html.UserData;v=s.viewState;if isnumeric(name)&&isequal(name,0),name='';end;name=char(name);
        regions=s.payload.data.evidenceLedger.regions;idx=find(strcmp({regions.name},name),1);
        assert(isempty(name)||~isempty(idx),'RocketAeroHeat:EvidenceView:Region','未知区域。');
        v.regionName=name;
        if ~isempty(idx)&&~isempty(regions(idx).faceIDs)
            v.selectedFace=regions(idx).faceIDs(1);v.position=mean(s.payload.data.vertices(s.payload.data.faces(v.selectedFace+1,:)+1,:),1);v.sectionX=v.position(1);
        end
        publish(v);
    end
    function inspectFace(i,position)
        s=html.UserData;d=s.payload.data;validateattributes(i,{'numeric'},{'scalar','integer','>=',1,'<=',size(d.faces,1)});
        if nargin<2||isempty(position),position=mean(d.vertices(d.faces(i,:)+1,:),1);end
        validateattributes(position,{'numeric'},{'real','finite','numel',3});v=s.viewState;v.selectedFace=i-1;v.position=reshape(position,1,3);v.sectionX=position(1);v.regionName=char(string(d.region(i)));publish(v);
    end
    function setRange(low,high)
        validateattributes([low high],{'numeric'},{'real','finite','numel',2});assert(low<high,'色标下限必须小于上限。');s=html.UserData;v=s.viewState;v.ranges.(v.fieldKey)=[low high];publish(v);
    end
    function setSectionX(x)
        s=html.UserData;validateattributes(x,{'numeric'},{'real','finite','scalar','>=',min(s.payload.data.vertices(:,1)),'<=',max(s.payload.data.vertices(:,1))});v=s.viewState;v.sectionX=x;publish(v);
    end
    function value=getState()
        s=html.UserData;value=s.viewState;value.revision=s.revision;value.selectedFaceID=s.viewState.selectedFace+1;
        value.rendered=s.rendered;value.ledgerIdentity=s.payload.data.evidenceLedger.identity;
    end
    function value=getData(),s=html.UserData;value=s.payload.data;end
    function value=getSection(),s=html.UserData;value=s.sectionInfo;end
    function receipt=exportTo(path)
        % 模块三：独立CSV只导出当前变量的全部原面元；选区不会筛掉未知区域。
        % SHA及全部变量的面积汇总写入小型旁车JSON，不在每一行重复。
        s=html.UserData;d=s.payload.data;ledger=d.evidenceLedger;key=s.viewState.fieldKey;f=ledger.fields(strcmp({ledger.fields.key},key));
        path=char(path);[folder,name,ext]=fileparts(path);if isempty(ext),ext='.csv';end
        assert(strcmpi(ext,'.csv'),'证据导出文件须使用 .csv 扩展名。');if isempty(folder),folder=pwd;end
        path=fullfile(folder,[name ext]);
        assert(isfolder(folder),'导出目录不存在。');summaryPath=fullfile(folder,[name '.json']);
        fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'无法写入证据CSV。');guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',char(65279));
        headers={'face_id','region','area_m2','field_key','value','unit','availability_code','availability','qualification_code','qualification','method','scope','reason','evidence_ids'};
        writeRow(fid,headers);availability={'未知','工程计算','研究计算','给定条件','超出适用域'};qualification={'未验证','验证失败','验证通过'};
        for i=1:size(d.faces,1)
            value='';if isfinite(f.values(i)),value=sprintf('%.17g',f.values(i));end
            scope=f.scope;if isfield(f,'scopeNames')&&~isempty(f.scopeNames),scope=cellText(f.scopeNames,f.scopeIndex(i)+1);end
            evidence=cellValue(f.evidenceIDs,i);if isempty(evidence),evidence='';else,evidence=jsonencode(evidence);end
            writeRow(fid,{sprintf('%d',i),cellText(d.region,i),sprintf('%.17g',d.faceArea(i)),key,value,f.unit, ...
                sprintf('%d',f.availability(i)),availability{f.availability(i)+1},sprintf('%d',f.qualification(i)),qualification{f.qualification(i)+1}, ...
                cellText(f.methodNames,f.methodIndex(i)+1),scope,cellText(f.reason,i),evidence});
        end
        clear guard;
        compact=compactLedger(ledger);v=s.viewState;
        summary=struct('schemaVersion',1,'kind','aerothermal-evidence-export','csvFile',[name ext], ...
            'csvSHA256',sha256File(path),'rows',size(d.faces,1),'fieldKey',key, ...
            'identity',ledger.identity,'evidenceLedger',compact,'viewState',v, ...
            'exportScope','当前变量的全部原始面元；选区不筛选CSV，原值不按资格裁剪');
        fid=fopen(summaryPath,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(summary));clear guard;
        receipt=struct('csvPath',path,'summaryPath',summaryPath,'rows',size(d.faces,1),'fieldKey',key,'csvSHA256',summary.csvSHA256);
    end
    function exportSectionTo(path)
        % 模块四：截面CSV来自当前网页的真实面/平面交线，不能沿用旧版本回执。
        s=html.UserData;assert(isstruct(s.sectionInfo)&&isfield(s.sectionInfo,'csv')&& ...
            strcmp(s.sectionInfo.fieldKey,s.viewState.fieldKey),'当前截面尚未完成绘制，请稍后导出。');
        fid=fopen(char(path),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',s.sectionInfo.csv);
    end
end

function range=initialRange(field)
% 模块五：初始范围只由全模型原数值确定；选择区域不重新归一化。
v=field.values(isfinite(field.values));if isempty(v),range=[0 1];return;end
range=[min(v),max(v)];if range(1)==range(2),delta=max(abs(range(1))*.01,1e-6);range=range+[-delta delta];end
end
function writeRow(fid,values)
quoted=cellfun(@(x)['"' strrep(char(string(x)),'"','""') '"'],values,'UniformOutput',false);fprintf(fid,'%s\r\n',strjoin(quoted,','));
end
function value=cellValue(values,i)
if iscell(values),value=values{i};else,value=values(i);end
end
function value=cellText(values,i)
value=char(string(cellValue(values,i)));
end
function value=compactLedger(ledger)
% 模块六：旁车保留身份/汇总/证据；逐面原值已经在独立CSV中，不重复大数组。
value=struct('schemaVersion',ledger.schemaVersion,'kind',ledger.kind,'identity',ledger.identity, ...
    'fields',struct([]),'regions',struct([]),'joint',withoutFaceIDs(ledger.joint), ...
    'acceptedEvidence',ledger.acceptedEvidence,'rejectedEvidence',ledger.rejectedEvidence);
for i=1:numel(ledger.fields)
    f=ledger.fields(i);entry=struct('key',f.key,'label',f.label,'unit',f.unit,'valueSHA256',f.valueSHA256,'areaSummary',f.areaSummary);
    if isempty(value.fields),value.fields=entry;else,value.fields(end+1)=entry;end
end
for i=1:numel(ledger.regions)
    r=ledger.regions(i);entry=struct('name',r.name,'area',r.area,'faceCount',numel(r.faceIDs),'fields',r.fields,'joint',withoutFaceIDs(r.joint));
    if isempty(value.regions),value.regions=entry;else,value.regions(end+1)=entry;end
end
end
function value=withoutFaceIDs(value)
names=fieldnames(value);remove=names(endsWith(names,'FaceIDs'));if ~isempty(remove),value=rmfield(value,remove);end
end
function hash=sha256File(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));digest=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),bytes=fread(fid,1024*1024,'*uint8');digest.update(typecast(bytes,'int8'));end
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
