function check_pair(out)
% 模块一：按原版本绑定两个保存场的几何、状态和源身份，进入工作台的实际区域页。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
engineeringFile=fullfile(root,'outputs','v104_validation','baseline02','current_case.mat');
researchFile=fullfile(root,'outputs','v167_validation','main01','result.mat');
d=load(engineeringFile,'current');current=d.current;
d=load(researchFile,'r');r=d.r;
assert(isequal(current.V,r.V)&&isequal(current.F,r.F)&&~r.productionQualified);
current.evidenceIdentity=struct('origin','engineering','sourceTraceable',true,'sourceFrame',struct( ...
    'kind','saved-engineering-surface','path',engineeringFile,'sha256',cm011FileSHA256(engineeringFile)));
r.geometryRegion=current.localGeometry.region;
r.geometryRegionSource='与研究V/F逐项相同的工程原面几何标签；没有转移工程流动资格';
fig=uifigure('Visible','off');clean=onCleanup(@()delete(fig)); %#ok<NASGU>
engineering=buildAerothermalViewerData(current);
view=createAerothermalEvidenceView(fig,engineering);
view.inspectFace(48,[]);
e=view.exportWallPairTo(fullfile(out,'engineering.csv'));
assert(e.rows==19752&&abs(e.joint.pressureHeatFluxEngineeringArea-.8467538430106993)<1e-11);
view.update(buildSU2EvidenceViewerData(r));
rs=view.exportWallPairTo(fullfile(out,'research.csv'));
assert(rs.rows==19752&&rs.joint.pressureHeatFluxEngineeringArea==0&&rs.joint.pressureHeatFluxQualifiedArea==0);
assert(strcmp(e.identity.geometrySHA256,rs.identity.geometrySHA256)&&~strcmp(e.identity.sourceSHA256,rs.identity.sourceSHA256));

% 模块二：从正常数据构建器生成小型快照，核对零值、负值、给定、缺项和历史保护。
original=view.getData();f=original.evidenceLedger.fields;
original.evidenceLedger.fields=f(~strcmp({f.key},'pressure'));
expectError(original,'RAH:WallPair:Fields',out);
original.evidenceLedger.fields=f;
original.evidenceLedger.fields(strcmp({f.key},'heatFlux')).history=[1 2];
expectError(original,'RAH:WallPair:History',out);
sample=engineering;
sample.evidenceIdentity=struct('origin','engineering','synthetic',true);
sample.resultCondition=struct('caption','三面软件测试');
sample.faces=sample.faces(1:3,:);sample.faceArea=sample.faceArea(1:3);sample.region=sample.region(1:3);
fields=sample.fields(ismember({sample.fields.key},{'pressure','heatFlux'}));
for j=1:numel(fields)
    fields(j).values=[0;-2;NaN];fields(j).status=[1;2;0];
    fields(j).methodNames={'方法甲','给定条件','未知'};fields(j).methodIndex=[0;1;2];
    fields(j).scopeNames={};fields(j).scopeIndex=[];
    fields(j).missingReason={'','','源记录缺项'};fields(j).outOfDomainFace=[false;false;true];
end
sample.fields=fields;sample.evidenceLedger=buildAerothermalEvidenceLedger(sample);
exportWallPairEvidence(sample,fullfile(out,'signed.csv'));
% 复用已有的六面合成证据，检验已验证、已计算与超域的独立判定。
fixture=fullfile(root,'outputs','v176_validation','unit01','fixture.json');
assert(strcmpi(cm011FileSHA256(fixture),'acbc9da78354199669a13e979230b9422b46856c443e2568926f1e533eac3f8f'));
proof=jsondecode(fileread(fixture));assert(proof.syntheticOnly);
small=proof.data;small.evidenceLedger=proof.qualifiedLedger;
exportWallPairEvidence(small,fullfile(out,'qualified_synthetic.csv'));
fid=fopen(fullfile(out,'view.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(struct('engineering',e.identity,'research',rs.identity,'selectedFace',48, ...
    'engineeringRows',e.rows,'researchRows',rs.rows,'missingFieldRejected',true,'historyRejected',true)));
end

function expectError(data,id,out)
path=fullfile(out,'rejected.csv');
try
    exportWallPairEvidence(data,path);error('测试输入应当被拒绝');
catch ex
    assert(strcmp(ex.identifier,id),ex.message);assert(~isfile(path));
end
end
