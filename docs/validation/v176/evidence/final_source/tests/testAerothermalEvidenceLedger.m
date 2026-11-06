function report=testAerothermalEvidenceLedger(out)
%TESTAEROTHERMAL EVIDENCE LEDGER 检验逐变量资格，绝不将合成证据用于实际物理认证。
% 模块一：六个真正不等面积三角面，计算值/给定值/未知/超域分别可辨。
% 模块二：真实读取有SHA的合成proof，检查交集、失败优先、重复及错误来源。
% 模块三：保存原始canonical输入和输出，供Python独立计算面积与集合。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(~isfolder(out),'测试证据目录已存在，拒绝覆盖');mkdir(out);
timer=tic; data=fixture(); before=data;
initial=buildAerothermalEvidenceLedger(data);
assert(isequaln(data,before),'台账修改了原始场量');
assert(initial.totalArea==63);
p=field(initial,'pressure');q=field(initial,'heatFlux');
assert(p.areaSummary.engineeringArea==39 && p.areaSummary.prescribedArea==8 && p.areaSummary.unknownArea==16);
assert(q.areaSummary.engineeringArea==15 && q.areaSummary.unknownArea==16 && q.areaSummary.outOfDomainArea==32);
assert(q.values(1)==0 && q.values(2)==-2 && isnan(q.values(5)) && q.values(6)==6);
assert(all(p.qualification==0) && q.qualification(6)==1);
assert(initial.joint.pressureHeatFluxQualifiedArea==0);

% 合成数值证据只是单元测试；所有记录都绑定完整几何、工况、来源与场值。
pPass=proofRecord(out,initial,'pressure',[0 1 2],'synthetic-pressure-pass',false);
qPass=proofRecord(out,initial,'heatFlux',[0 1],'synthetic-heat-zero-negative-pass',false);
qConflictPass=proofRecord(out,initial,'heatFlux',2,'synthetic-heat-face2-pass',false);
qFail=proofRecord(out,initial,'heatFlux',2,'synthetic-heat-face2-fail',true);
records=[pPass qPass qConflictPass qFail];
qualified=buildAerothermalEvidenceLedger(data,records);
p=field(qualified,'pressure');q=field(qualified,'heatFlux');
assert(p.areaSummary.qualifiedArea==7 && q.areaSummary.qualifiedArea==3);
assert(q.qualification(3)==1,'不同证据冲突时失败必须优先');
assert(qualified.joint.pressureHeatFluxQualifiedArea==3);
assert(qualified.joint.pressureHeatFluxQualifiedArea~=mean([p.areaSummary.qualifiedArea q.areaSummary.qualifiedArea]));
assert(isequaln(data,before) && isequaln(p.values,data.fields(1).values));
assert(isequaln(q.values,data.fields(2).values),'超域有限值不能被补零或静默裁掉');
implicitData=data;implicitData.evidenceRecords=records;
implicit=buildAerothermalEvidenceLedger(implicitData);
assert(isequaln(implicit.joint,qualified.joint) && implicit.acceptedEvidenceCount==4);
assert(isequaln(field(implicit,'heatFlux').qualification,q.qualification));
nanOutsideData=data;nanOutsideData.fields(2).outOfDomainFace(5)=true;
nanOutside=buildAerothermalEvidenceLedger(nanOutsideData);item=field(nanOutside,'heatFlux');
assert(isnan(item.values(5)) && item.status(5)==0 && item.availability(5)==4 && item.qualification(5)==1);
assert(item.areaSummary.outOfDomainArea==48 && item.areaSummary.unknownArea==0);

duplicate=buildAerothermalEvidenceLedger(data,[records pPass]);
assert(duplicate.duplicateEvidenceCount==1 && duplicate.acceptedEvidenceCount==4);
assert(field(duplicate,'pressure').areaSummary.qualifiedArea==7);
assert(duplicate.joint.pressureHeatFluxQualifiedArea==3);

% 相同id但内容冲突，整组拒绝，不能先接纳一条再留下资格。
conflict=proofRecord(out,initial,'pressure',[0 1 2],'synthetic-pressure-pass',true,'same_id_conflict');
conflicted=buildAerothermalEvidenceLedger(data,[pPass conflict]);
assert(conflicted.rejectedEvidenceCount==2 && field(conflicted,'pressure').areaSummary.qualifiedArea==0);
rejections=struct('name',{},'count',{});
rejections(end+1)=struct('name','same-id-conflict-rejects-whole-group','count',2);

% 每个错误同时修改record与proof，避免仅由二者不一致挡住真实绑定反例。
for name={'geometrySHA256','conditionSHA256','sourceSHA256','valueSHA256'}
    changed=pPass;changed.bindings.(name{1})=repmat('f',1,64);
    changed=rewriteProof(out,changed,['wrong_' name{1}]);
    bad=buildAerothermalEvidenceLedger(data,changed);
    assert(bad.rejectedEvidenceCount==1 && field(bad,'pressure').areaSummary.qualifiedArea==0);
    rejections(end+1)=struct('name',['wrong-' name{1}],'count',1); %#ok<AGROW>
end
for mode=1:3
    changed=pPass;
    if mode==1,changed.faceIDs=[1 0];label='face-order';
    elseif mode==2,changed.faceIDs=[0 0];label='duplicate-face';
    else,changed.faceIDs=6;label='outside-original-face';end
    changed=rewriteProof(out,changed,label);bad=buildAerothermalEvidenceLedger(data,changed);
    assert(bad.rejectedEvidenceCount==1 && field(bad,'pressure').areaSummary.qualifiedArea==0);
    rejections(end+1)=struct('name',label,'count',1); %#ok<AGROW>
end
changed=pPass;changed.proofSHA256=repmat('0',1,64);
bad=buildAerothermalEvidenceLedger(data,changed);
assert(bad.rejectedEvidenceCount==1);rejections(end+1)=struct('name','wrong-proof-file-sha','count',1);

% 泛型passed或缺失热边界/近壁规则不能替代逐变量数值证据。
changed=qPass;proof=jsondecode(fileread(changed.proofFile));proof=rmfield(proof,'criteria');proof.passed=true;
changed=saveChangedProof(out,changed,proof,'generic_passed_without_criteria');
bad=buildAerothermalEvidenceLedger(data,changed);
assert(bad.rejectedEvidenceCount==1);rejections(end+1)=struct('name','generic-passed-no-criteria','count',1);
changed=qPass;proof=jsondecode(fileread(changed.proofFile));
proof.criteria=proof.criteria(~strcmp({proof.criteria.name},'near-wall-resolution'));
changed=saveChangedProof(out,changed,proof,'heat_missing_nearwall');
bad=buildAerothermalEvidenceLedger(data,changed);
assert(bad.rejectedEvidenceCount==1);rejections(end+1)=struct('name','heat-without-nearwall-rule','count',1);
notSynthetic=data;notSynthetic.evidenceIdentity.synthetic=false;
bad=buildAerothermalEvidenceLedger(notSynthetic,pPass);
assert(bad.rejectedEvidenceCount==1);rejections(end+1)=struct('name','synthetic-proof-on-real-data','count',1);

changedData=data;changedData.resultCondition.Mach=6;
bad=buildAerothermalEvidenceLedger(changedData,pPass);
assert(bad.rejectedEvidenceCount==1);rejections(end+1)=struct('name','stale-condition-data','count',1);
changedData=data;changedData.evidenceIdentity.sourceFrame.caseID='different-synthetic-source';
bad=buildAerothermalEvidenceLedger(changedData,pPass);
assert(bad.rejectedEvidenceCount==1);rejections(end+1)=struct('name','stale-source-data','count',1);
% 原值不变但算法、状态或逐变量来源改变，旧proof也必须失效。
changedData=data;changedData.fields(1).methodNames={'不同的合成方法'};
bad=buildAerothermalEvidenceLedger(changedData,pPass);
assert(bad.rejectedEvidenceCount==1 && ~strcmp(field(bad,'pressure').valueSHA256,field(initial,'pressure').valueSHA256));
rejections(end+1)=struct('name','same-values-changed-method','count',1);
changedData=data;changedData.fields(1).status(1)=2;
bad=buildAerothermalEvidenceLedger(changedData,pPass);
assert(bad.rejectedEvidenceCount==1 && ~strcmp(field(bad,'pressure').valueSHA256,field(initial,'pressure').valueSHA256));
rejections(end+1)=struct('name','same-values-changed-status','count',1);
changedData=data;changedData.fields(1).evidenceOrigin='research';
bad=buildAerothermalEvidenceLedger(changedData,pPass);
assert(bad.rejectedEvidenceCount==1 && ~strcmp(field(bad,'pressure').valueSHA256,field(initial,'pressure').valueSHA256));
rejections(end+1)=struct('name','same-values-changed-field-origin','count',1);
changedData=data;changedData.faces([1 2],:)=changedData.faces([2 1],:);
changedData.faceArea([1 2])=changedData.faceArea([2 1]);
bad=buildAerothermalEvidenceLedger(changedData,pPass);
assert(bad.rejectedEvidenceCount==1);rejections(end+1)=struct('name','original-face-order-changed','count',1);

% 给定Tw和未知固体Tw不能被proof升级；研究数据也不能计入工程面积。
wallProof=proofRecord(out,initial,'wallInput',[0 1 2 3 4 5],'synthetic-wall-given-no-upgrade',false);
unknownProof=proofRecord(out,initial,'wallSolved',[0 1 2 3 4 5],'synthetic-wall-unknown-no-upgrade',false);
protected=buildAerothermalEvidenceLedger(data,[wallProof unknownProof]);
assert(field(protected,'wallInput').areaSummary.qualifiedArea==0);
assert(field(protected,'wallInput').areaSummary.prescribedArea==63);
assert(field(protected,'wallSolved').areaSummary.qualifiedArea==0);
assert(field(protected,'wallSolved').areaSummary.unknownArea==63);
researchData=data;researchData.evidenceIdentity.origin='research';
research=buildAerothermalEvidenceLedger(researchData);
assert(field(research,'pressure').areaSummary.engineeringArea==0);
assert(field(research,'pressure').areaSummary.researchArea==39);
assert(field(research,'heatFlux').areaSummary.researchArea==15);
assert(research.joint.pressureHeatFluxQualifiedArea==0);
legacyData=rmfield(data,'evidenceIdentity');legacy=buildAerothermalEvidenceLedger(legacyData);
assert(field(legacy,'pressure').areaSummary.qualifiedArea==0);
% 历史SU2的均匀壁温可为标量；仅该边界输入允许展开，原物理场不许广播。
surface=struct('V',data.vertices(1:6,:),'F',data.faces(1:2,:)+1, ...
    'pressureRaw',[10;20],'wallTemperatureInput',300);
adapted=buildSU2EvidenceViewerData(surface);wall=field(adapted.evidenceLedger,'wallInput');
assert(isequal(wall.values,[300;300]) && isequal(wall.status,[2;2]));
surface.wallTemperatureInput=NaN;adapted=buildSU2EvidenceViewerData(surface);
wall=field(adapted.evidenceLedger,'wallInput');assert(all(isnan(wall.values)) && all(wall.status==0));
surface.pressureRaw=10;rejected=false;
try,buildSU2EvidenceViewerData(surface);catch ex,rejected=strcmp(ex.identifier,'RocketAeroHeat:SU2Evidence:Contract');end
assert(rejected,'原面压力场的标量不得被偷偷广播');

payload=struct('schema','RAH.EvidenceLedger.IndependentFixture.v1','syntheticOnly',true, ...
    'newQualifiedPhysicalArea',0,'data',data,'initialLedger',initial,'records',records, ...
    'qualifiedLedger',qualified,'duplicateLedger',duplicate,'conflictedLedger',conflicted, ...
    'implicitData',implicitData,'implicitLedger',implicit,'nanOutsideData',nanOutsideData,'nanOutsideLedger',nanOutside, ...
    'protectedLedger',protected,'researchData',researchData,'researchLedger',research, ...
    'legacyData',legacyData,'legacyLedger',legacy,'rejections',rejections);
writeJSON(fullfile(out,'fixture.json'),payload);
report=struct('passed',true,'syntheticOnly',true,'faces',6,'knownExactAreas',[1 2 4 8 16 32], ...
    'totalArea',63,'pressureQualifiedArea',7,'heatQualifiedArea',3,'jointQualifiedArea',3, ...
    'rejectionCases',numel(rejections),'rejectedEvidenceRecords',sum([rejections.count]), ...
    'zeroAndNegativeHeatPreserved',true,'unknownAndOutOfDomainSeparated',true, ...
    'researchNotCountedAsEngineering',true,'rawDataUnchanged',isequaln(data,before), ...
    'legacyScalarWallInputCases',3, ...
    'newPhysicalArea',0,'newCFDTasks',0,'elapsedSeconds',toc(timer));
writeJSON(fullfile(out,'report.json'),report);
fprintf('Evidence ledger fixture: 6 unequal faces, p=7/q=3/intersection=3; %d rejection cases; no physical qualification.\n',numel(rejections));
end

function data=fixture()
areas=[1;2;4;8;16;32];vertices=zeros(18,3);faces=zeros(6,3);
for j=1:6
    ids=3*j-2:3*j;vertices(ids,:)=[10*j 0 0;10*j+2*areas(j) 0 0;10*j 1 0];faces(j,:)=ids-1;
end
region={'nose';'sidewall';'sidewall';'base';'base';'appendage'};
data=struct('kind','aerothermal-model','schemaVersion',1,'vertices',vertices,'faces',faces, ...
    'faceArea',areas,'region',{region},'fields',struct([]), ...
    'resultCondition',struct('Mach',5,'alpha',20,'beta',0,'altitude',30000,'Twall',300,'caption','明确合成工况'), ...
    'evidenceIdentity',struct('origin','engineering','sourceTraceable',true,'synthetic',true, ...
        'sourceFrame',struct('caseID','synthetic-six-unequal-triangles','frame',0)));
data.fields=[makeField('pressure','压力','kPa',[10;20;30;40;NaN;60],[1;1;1;2;0;1]), ...
    makeField('heatFlux','热流','kW/m^2',[0;-2;3;4;NaN;6],[1;1;1;1;0;1]), ...
    makeField('wallInput','给定壁温','K',300*ones(6,1),2*ones(6,1)), ...
    makeField('wallSolved','固体壁温解','K',nan(6,1),zeros(6,1))];
data.fields(2).outOfDomainFace(6)=true;
end

function value=makeField(key,label,unit,values,status)
value=struct('key',key,'label',label,'unit',unit,'values',values,'status',status, ...
    'methodNames',{{['合成测试方法-' key]}},'methodIndex',zeros(6,1), ...
    'scope','仅验证台账合同，不是物理方法认证','outOfDomainFace',false(6,1));
end

function value=field(ledger,key)
value=ledger.fields(strcmp({ledger.fields.key},key));assert(isscalar(value));
end

function record=proofRecord(out,ledger,variable,faceIDs,id,fail,label)
if nargin<7,label=id;end
item=field(ledger,variable);bindings=ledger.identity;
bindings=struct('geometrySHA256',bindings.geometrySHA256,'conditionSHA256',bindings.conditionSHA256, ...
    'sourceSHA256',bindings.sourceSHA256,'valueSHA256',item.valueSHA256);
names={'condition-match','variable-reference-comparison','numerical-accuracy','model-region-applicability'};
if strcmp(variable,'heatFlux'),names=[names {'thermal-boundary-properties','near-wall-resolution'}];end
criteria=repmat(struct('name','','metric','synthetic-independent-identity','observed',0, ...
    'operator','==','threshold',0,'unit','1'),1,numel(names));
for j=1:numel(names),criteria(j).name=names{j};end
if fail,criteria(2).observed=1;end
proof=struct('schema','RAH.AerothermalQualificationProof.v1','evidenceID',id,'variable',variable, ...
    'faceIDs',faceIDs,'bindings',bindings,'basis','synthetic-fixture','applicability','仅该合成六面元输入', ...
    'reference',struct('type','synthetic-fixture','uri','synthetic-fixture://unequal-area-unit-test', ...
        'sha256',repmat('a',1,64)),'criteria',criteria);
path=fullfile(out,[label '.json']);writeJSON(path,proof);
record=struct('schema','RAH.AerothermalEvidenceRecord.v1','id',id,'variable',variable,'faceIDs',faceIDs, ...
    'bindings',bindings,'proofFile',absolute(path),'proofSHA256',fileSHA(path));
end

function record=rewriteProof(out,record,label)
proof=jsondecode(fileread(record.proofFile));proof.bindings=record.bindings;proof.faceIDs=record.faceIDs;
record=saveChangedProof(out,record,proof,label);
end

function record=saveChangedProof(out,record,proof,label)
path=fullfile(out,[label '.json']);writeJSON(path,proof);
record.proofFile=absolute(path);record.proofSHA256=fileSHA(path);
end

function value=absolute(path)
value=char(java.io.File(path).getCanonicalPath());
end

function value=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function writeJSON(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'ConvertInfAndNaN',true));
end
