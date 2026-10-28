function testAerothermalViewerTemperatures(out,syntheticOnly)
%TESTAEROTHERMALVIEWERTEMPERATURES 校核单工况温度场的值、来源与未解保护。
% 可选第二参数仅运行合成查看器合同，避免为界面元数据验证重复求解物理。
if nargin<2,syntheticOnly=false;end
% 模块一：受控六面元软件样本分别覆盖层流、湍流、球鼻、脱体闭合
% 和未解；该样本检验数据合同，不作为任何气动热实验验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
m=struct('V',zeros(18,3),'F',reshape(1:18,3,6)');
for i=1:6
    j=(i-1)*3+(1:3);m.V(j,:)=[i 0 0;i 2*i 0;i 0 1];
end
ids=uint8([1;2;3;4;5;0]);Pr=.72;Te=400;U=1000;
Tlam=airRecoveryTemperature(Te,U,sqrt(Pr));
Tturb=airRecoveryTemperature(Te,U,Pr^(1/3));
Tnose=airRecoveryTemperature(220,1400,1);
Taw=[Tlam;Tturb;Tnose;Tlam;Tturb;999];
m.heating=struct('heatFlux',(1:6)'*1000,'eligibleFace',ids>0, ...
    'faceMethodID',ids,'faceMethodNames',{{'未计算','热流层流','热流湍流', ...
    '球鼻热流','脱体层流热流','脱体湍流热流'}}, ...
    'validity','仅软件来源样本','recoveryTemperature',Taw, ...
    'wallTemperatureByFace',[300;320;340;360;NaN;400], ...
    'Prandtl',Pr,'thermalPropertyMethod','thermallyPerfect', ...
    'thermalPropertyScope','干空气焓模型 80~2500 K，不含解离', ...
    'boundaryLayerMethod','Thwaites 热流积分', ...
    'turbulentMethod','Head 热流积分','stagnationMethod','Sutton-Graves 热流', ...
    'noseMethod','驻点热流余弦延拓','noseScope','验证球形鼻端与入射角范围');
m.stagnationHeating=struct('adiabaticWallTemperature',Tnose, ...
    'thermalPropertyMode','thermallyPerfect', ...
    'shockStateMethod','法向激波总焓守恒（样本）');
m.boundaryLayerEdge=struct('pressure',ones(6,1)*10000, ...
    'temperature',[400;400;500;400;400;999],'algebraicallyValid',ids>0, ...
    'faceMethodID',uint8([2;3;5;5;4;0]), ...
    'faceMethodNames',{{'无解','近切向','切锥','切楔','膨胀','前体闭合'}}, ...
    'thermalPropertyModeByFace',uint8([2;2;2;2;1;0]), ...
    'thermalPropertyModeNames',{{'无解','定比热','NASA7 干空气'}}, ...
    'thermalExpansionOutOfRangeFace',[false;false;false;false;true;false], ...
    'scope','局部外流适用性样本');
m.pressure=struct('model','表面压力样本','staticPressure',ones(6,1)*1234.5, ...
    'Cp',zeros(6,1),'isPrescribedLeeward',[true;false;false;false;false;false], ...
    'freestream',struct('Mach',5,'alpha',20,'beta',0,'altitude',30000, ...
    'pressure',1234.5,'temperature',226.5,'unitFlowDirection',[-1 0 0]));
% 故意存入与结果不一致的编辑中配置；查看器只能使用冻结的结果自由流。
m.cfg=struct('flight',struct('Mach',99,'alpha',88,'altitude',86000,'Twall',2000));
d=buildAerothermalViewerData(m);r=get(d,'recoveryTemperature');
assert(isequaln(r.values,[Taw(1:5);NaN]) && isequal(r.status,[1;1;1;1;1;0]));
names=decode(r);scopes=decodeScopes(r);
assert(contains(names{1},'层流 r=√Pr') && contains(names{2},'湍流 r=Pr^(1/3)'));
assert(contains(names{1},'NASA7') && contains(names{3},'总温 T0'));
assert(contains(names{3},'法向激波总焓守恒') && contains(scopes{3},'不是各面元局部'));
assert(contains(scopes{4},'脱体') && ~contains(scopes{1},'脱体'));
assert(~any(contains(names,{'Sutton-Graves','Thwaites','Head','热流积分'})));
assert(contains(names{6},'未解') && abs(r.missingAreaFraction-6/21)<1e-12);
% 恢复温度独立检查焓平衡；查看器输出不能改动真实恢复温度。
[~,he]=airThermalProperties(Te);[~,hr]=airThermalProperties(r.values([1 2]));
assert(max(abs((hr-he)-[sqrt(Pr);Pr^(1/3)]*U^2/2))<1e-5);
w=get(d,'wallInput');assert(isequal(w.status,[2;2;2;2;0;2]) && isnan(w.values(5)));
solved=get(d,'wallSolved');assert(all(solved.status==0) && all(isnan(solved.values)));
assert(contains(solved.scope,'独立热平衡') && solved.computedAreaFraction==0);
edge=get(d,'edgeTemperature');edgeMethods=decode(edge);edgeScopes=decodeScopes(edge);
assert(contains(edgeMethods{1},'NASA7 干空气') && contains(edgeMethods{5},'定比热'));
assert(contains(edgeScopes{5},'流线占位') && isnan(edge.values(6)));
assert(edge.status(5)==0 && isnan(edge.values(5)), ...
    '域外定比热流向占位不能显示为已计算的外缘温度。');
edgePressure=get(d,'edgePressure');
assert(edgePressure.status(5)==0 && isnan(edgePressure.values(5)));
assert(any(contains(decode(edgePressure),'未解')));
assert(contains(edge.missingReason{5},'本面元外缘热力状态未知'));
assert(contains(edgeMethods{5},'未接受的尝试模型') && contains(edgeMethods{5},'定比热'));
assert(contains(w.missingReason{5},'输入壁温未知') && ~contains(w.missingReason{5},'上游'));
assert(strcmp(solved.availability,'notProvidedInSnapshot') && contains(solved.missingReason{1},'未执行壁温求解'));
assert(d.resultCondition.Mach==5 && d.resultCondition.altitudeM==30000);
assert(d.resultCondition.freestreamPressurePa==1234.5 && d.resultCondition.freestreamTemperatureK==226.5);
assert(strcmp(d.resultCondition.wallInputType,'perFace') && d.resultCondition.wallUnknownFaces==1);
assert(isnan(d.resultCondition.wallUniformValueK) && ~contains(d.caseNote,'2000 K'));
assert(contains(d.caseNote,'高度 30 km') && contains(d.caseNote,'p∞ 1234.5 Pa') && contains(d.caseNote,'不是壁温解'));
% 本面外缘未知与未知上游导致热边界缺失，局部外缘有效仍不能恢复热流。
guarded=m;guarded.heating.eligibleFace([5 6])=false;
guarded.heating.heatFlux([5 6])=NaN;guarded.heating.recoveryTemperature([5 6])=NaN;
guarded.boundaryLayerEdge.algebraicallyValid(6)=true;guarded.boundaryLayerEdge.temperature(6)=400;
guarded.continuumRegime=struct('directStateUnresolvedFace',[false;false;false;false;true;false], ...
    'downstreamStateUnresolvedFace',[false;false;false;false;false;true]);
guardedData=buildAerothermalViewerData(guarded);hg=get(guardedData,'heatFlux');eg=get(guardedData,'edgeTemperature');
assert(contains(hg.missingReason{5},'本面元外缘') && contains(hg.missingReason{6},'上游外缘'));
assert(contains(hg.missingReason{6},'即使本面外缘有效') && eg.status(6)==1 && hg.status(6)==0);
assert(isequaln(hg.values,guarded.heating.heatFlux/1000) && all(hg.values(1:4)==(1:4)'));
uniform=m;uniform.heating.wallTemperatureByFace(:)=300;
du=buildAerothermalViewerData(uniform);assert(strcmp(du.resultCondition.wallInputType,'uniform'));
assert(du.resultCondition.wallUniformValueK==300 && contains(du.caseNote,'统一壁温边界输入 300 K'));
noRecord=rmfield(uniform,'pressure');dn=buildAerothermalViewerData(noRecord);
assert(isnan(dn.resultCondition.altitudeM) && contains(dn.caseNote,'高度 未记录 km'));
assert(~contains(dn.caseNote,'86 km')); % 不从待编辑cfg推断缺失工况。
% 显式定比热与旧结果缺失来源的情况也不得误标成干空气焓恢复。
old=m;old.heating.thermalPropertyMethod='caloricallyPerfect';
c=get(buildAerothermalViewerData(old),'recoveryTemperature');
cn=decode(c);assert(all(contains(cn([1 2 4 5]),'cp=gamma*R/(gamma-1)')));
old.heating=rmfield(old.heating,{'thermalPropertyMethod','Prandtl'});
c=get(buildAerothermalViewerData(old),'recoveryTemperature');cn=decode(c);
assert(contains(cn{1},'热物性模式未记录') && contains(cn{1},'Pr 数值未记录'));
% 无局部编号的外缘仍使用原始整体方法，不丢失真实来源。
old.boundaryLayerEdge=rmfield(old.boundaryLayerEdge,{'faceMethodID','faceMethodNames'});
old.boundaryLayerEdge.method='原始全局外流方法';
c=get(buildAerothermalViewerData(old),'edgeTemperature');
assert(all(contains(decode(c),'原始全局外流方法')));
% 模块二：经典锥柱真实内部外流/热流链路；合成模式明确跳过这项物理复算。
mesh=[];cfg=[];realData=[];
if ~syntheticOnly
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=5;
cfg.flight.beta=0;cfg.flight.Twall=300;
mesh=buildCanonicalBody('coneCylinder',struct('resolution',24));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);mesh=calcFlightStationState(mesh,cfg.flight,cfg);
realData=buildAerothermalViewerData(mesh);r=get(realData,'recoveryTemperature');
expected=mesh.heating.recoveryTemperature;expected(~mesh.heating.eligibleFace)=NaN;
assert(isequaln(r.values,expected));
names=decode(r);mask=ismember(double(mesh.heating.faceMethodID),[1 4]);
assert(any(mask) && all(contains(names(mask),'层流 r=√Pr')));
mask=ismember(double(mesh.heating.faceMethodID),[2 5]);
if any(mask),assert(all(contains(names(mask),'湍流 r=Pr^(1/3)')));end
assert(~any(contains(names,{'Sutton-Graves','Thwaites','Head'})));
solved=get(realData,'wallSolved');assert(all(isnan(solved.values)));
input=get(realData,'wallInput');assert(all(input.status==2) && all(input.values==300));
end
if nargin>0
    if ~exist(out,'dir'),mkdir(out);end
    save(fullfile(out,'single_time_temperatures.mat'),'mesh','cfg','realData');
    file=fopen(fullfile(out,'temperature_provenance_fixture.json'),'w','n','UTF-8');
    clean=onCleanup(@()fclose(file)); %#ok<NASGU>
    fprintf(file,'%s',jsonencode(d,'ConvertInfAndNaN',true));
    clear clean;
    file=fopen(fullfile(out,'viewer_field_state_fixture.json'),'w','n','UTF-8');
    assert(file>=0,'不能写入逐变量状态样本');clean=onCleanup(@()fclose(file)); %#ok<NASGU>
    fprintf(file,'%s',jsonencode(guardedData,'ConvertInfAndNaN',true));
end
if syntheticOnly
    fprintf('Synthetic viewer contract passed: six unequal-area faces; per-field missing/attempted provenance; frozen freestream and actual wall-input identity. No physical pipeline recomputed.\n');
else
fprintf(['Single-time temperature provenance passed: 6-method unequal-area fixture; ' ...
    'NASA7 enthalpy residual <1e-5 J/kg; %d-face real cone/cylinder chain; ' ...
    'input wall temperature prescribed, absent wall solution NaN.\n'],size(mesh.F,1));
end
end
function field=get(data,key)
field=data.fields(strcmp({data.fields.key},key));assert(isscalar(field));
end
function names=decode(field)
names=reshape(field.methodNames(field.methodIndex+1),[],1);
end
function scopes=decodeScopes(field)
scopes=reshape(field.scopeNames(field.scopeIndex+1),[],1);
end
