function report=testTmrAxisymmetricSutherlandView(folder,evidenceFile)
%TESTTMRAXISYMMETRICSUTHERLANDVIEW 核对真实变输运稳态页面及旧常数页回放。
% 模块一：读取实际Sutherland求稳态证据，重用几何/场量/曲面账及原18项拒绝。
% 模块二：参考muRef/Tref/S、给定Tw、实际局部墙mu/k与未知Te/Tr各自保留身份。
% 模块三：旧常数JSON除当前源码SHA和证据路径外逐值相同，HTML及方法文案原样保留。
% 模块四：新增输运元数据/局部系数篡改须拒绝；不运行长段求解或修改共享UI。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'tests'),fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v64_validation','tmr_sutherland_view');end
if nargin<2
    evidenceFile=pick({fullfile(root,'outputs','v64_validation','tmr_sutherland','tmr_sutherland_evidence.mat'), ...
        fullfile(root,'validation','v64','tmr_sutherland','tmr_sutherland_evidence.mat')});
end
assert(exist(evidenceFile,'file')==2,'AxisFV:SutherlandViewEvidence','需要真实Sutherland稳态证据。');
sourceSHA=fileSHA(evidenceFile);saved=load(evidenceFile,'mesh','implicit','gas','boundary','report');
assert(saved.report.steadyConverged&&saved.implicit.steadyConverged&& ...
    isfield(saved.gas,'transport')&&strcmp(saved.gas.transport.model,'sutherland'), ...
    'AxisFV:SutherlandViewEvidence','未收敛或常输运证据不能充当变输运稳态。');
uiPaths={'web/index.html','web/viewer.css','web/viewer.js','web/section.js','web/export.js'};
uiBefore=cellfun(@(p)fileSHA(fullfile(root,p)),uiPaths,'UniformOutput',false);

%% 一、原18项拒绝与几何/三角曲面账直接在真实变输运结果上重新执行
directed=testTmrAxisymmetricSteadyView(folder,evidenceFile);
assert(numel(directed.refusals)==18,'AxisFV:SutherlandViewRefusals','原18项拒绝必须全部保留。');
data=jsondecode(fileread(directed.view.jsonFile));mesh=saved.mesh;result=saved.implicit;
gas=axisymmetricFVGas(saved.gas);w=result.finalInfo.wall;map=data.wallBinding.sourceWallIndex;
tr=gas.transport;
assert(tr.referenceViscosity==1.716e-5&&tr.referenceTemperature==491.6*(5/9)&& ...
    tr.sutherlandTemperature==198.6*(5/9),'AxisFV:SutherlandViewReference','须保留官方Rankine参考温度的SI换算。');
assert(isfield(data,'transport')&&data.researchIdentity.localVariableTransport&& ...
    contains(data.researchIdentity.kind,'Sutherland')&&~contains(data.researchIdentity.kind,'constant-property'));
transport=data.transport;
assert(isequaln(transport.modelInput,tr)&&transport.referenceViscosityPaS==gas.mu&& ...
    transport.referenceTemperatureK==tr.referenceTemperature&&transport.referenceConductivityWmK==gas.k&& ...
    transport.constantPrandtl==gas.Pr&&isfield(transport,'parameterProvenance'));
assert(transport.parameterProvenance.referenceTemperatureRankine==491.6&& ...
    transport.parameterProvenance.sutherlandTemperatureRankine==198.6);
assert(contains(transport.scalarGasMuMeaning,'not current')&&contains(transport.scalarGasKMeaning,'not current'));
% 独立按明确公式核对元数据映射，防止把输入muRef直接贴到所有墙面。
expectedMu=tr.referenceViscosity*(w.temperatureInput/tr.referenceTemperature).^1.5.* ...
    ((tr.referenceTemperature+tr.sutherlandTemperature)./(w.temperatureInput+tr.sutherlandTemperature));
expectedK=expectedMu*(gas.cp/gas.Pr);
assert(max(abs(w.viscosity-expectedMu)./expectedMu)<5e-14&& ...
    max(abs(w.thermalConductivity-expectedK)./expectedK)<5e-14&&all(abs(w.viscosity/gas.mu-1)>1e-4));
assert(isequaln(transport.wallViscosityPaS,w.viscosity)&& ...
    isequaln(transport.wallThermalConductivityWmK,w.thermalConductivity)&& ...
    isequaln(data.wallBinding.sourceViscosityPaS,w.viscosity(map))&& ...
    isequaln(data.wallBinding.sourceThermalConductivityWmK,w.thermalConductivity(map))&& ...
    isequaln(data.wallBinding.transportEvaluationTemperatureInputK,w.temperatureInput(map)));
assert(isequaln(data.resultCondition.transportInput,tr)&& ...
    isequaln(data.resultCondition.wallViscosityPaS,w.viscosity)&& ...
    isequaln(data.resultCondition.wallThermalConductivityWmK,w.thermalConductivity));
assert(contains(data.caseNote,'Sutherland局部变输运')&&contains(data.caseNote,'muRef=')&& ...
    contains(data.caseNote,'参考T=')&&contains(data.caseNote,'mu(Tw)')&&~contains(data.caseNote,'常物性'));
html=fileread(directed.view.htmlFile);
assert(contains(html,'Sutherland局部变输运')&&~contains(html,'常物性层流'));
keys={data.fields.key};assert(numel(keys)==6&&isequal(keys, ...
    {'heatFlux','pressure','wallInput','wallSolved','edgeTemperature','recoveryTemperature'}));
q=data.fields(strcmp(keys,'heatFlux'));p=data.fields(strcmp(keys,'pressure'));Tw=data.fields(strcmp(keys,'wallInput'));
assert(contains(q.methodNames{1},'k(Tw)')&&contains(q.methodNames{1},'Sutherland')&& ...
    contains(p.methodNames{1},'EOS压力迹')&&contains(Tw.methodNames{1},'边界输入'));
for j=1:numel(data.fields)
    assert(all(contains(data.fields(j).scopeNames,'mu(Tw)'))&&all(contains(data.fields(j).scopeNames,'muRef=')));
end
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    field=data.fields(strcmp(keys,key{1}));assert(all(isnan(field.values))&&all(field.status==0));
end
paths={data.researchIdentity.currentSourceBindings.path};
assert(ismember('regionalFlow/axisymmetric/axisymmetricFVTransport.m',paths)&& ...
    ismember('regionalFlow/axisymmetric/createTmrAxisymmetricAirGas.m',paths));

%% 二、旧常数结果完整重放；只剔除本就随当前源码/解包路径改变的来源字段
legacyEvidence=pick({fullfile(root,'outputs','v63_validation','tmr_hlle_continue1','tmr_implicit_budget.mat'), ...
    fullfile(root,'validation','v63','tmr_hlle_continue1','tmr_implicit_budget.mat')});
legacyJSON=pick({fullfile(root,'outputs','v63_validation','tmr_steady_view','web','results.json'), ...
    fullfile(root,'validation','v63','tmr_steady_view','web','results.json')});
legacySHA=fileSHA(legacyEvidence);legacyFrozen=jsondecode(fileread(legacyJSON));
constant=testTmrAxisymmetricSteadyView(fullfile(folder,'legacy_constant'),legacyEvidence);
constantData=jsondecode(fileread(constant.view.jsonFile));
assert(numel(constant.refusals)==18&&~isfield(constantData,'transport')&& ...
    ~isfield(constantData.wallBinding,'sourceViscosityPaS')&&~isfield(constantData.resultCondition,'transportInput'));
assert(isequaln(withoutCurrentOrigin(constantData),withoutCurrentOrigin(legacyFrozen)), ...
    'AxisFV:SutherlandViewConstantRegression','旧常数JSON发生了来源SHA/路径之外的改变。');
assert(strcmp(fileSHA(constant.view.htmlFile),fileSHA(fullfile(fileparts(legacyJSON),'index.html'))), ...
    'AxisFV:SutherlandViewConstantRegression','旧常数页面HTML/文案发生改变。');

%% 三、变输运额外拒绝；每个失败必须发生在创建页面目录之前
refusals=struct('label',{},'identifier',{});badFolder=[tempname(folder) '_sutherland_must_not_export'];
bad=result;bad.gas=rmfield(bad.gas,'transport');reject('removed-local-transport',bad);
bad=result;bad.gas.transport.referenceTemperature=bad.gas.transport.referenceTemperature+1;reject('reference-temperature',bad);
bad=result;bad.gas.transport.sutherlandTemperature=bad.gas.transport.sutherlandTemperature+1;reject('sutherland-temperature',bad);
bad=result;bad.gas.transport.referenceViscosity=1.01*bad.gas.transport.referenceViscosity;reject('reference-mu-inconsistent',bad);
bad=result;bad.gas.k=1.01*bad.gas.k;reject('reference-k-inconsistent',bad);
bad=result;bad.finalInfo.wall.viscosity(1)=1.01*bad.finalInfo.wall.viscosity(1);reject('wall-local-mu',bad);
bad=result;bad.finalInfo.wall.thermalConductivity(1)=1.01*bad.finalInfo.wall.thermalConductivity(1);reject('wall-local-k',bad);
bad=result;bad.finalInfo.transport.faceMu(w.faceID(1))=gas.mu;reject('reference-mu-as-face-mu',bad);
bad=result;bad.finalInfo.transport.faceK(w.faceID(1))=gas.k;reject('reference-k-as-face-k',bad);
bad=result;bad.finalInfo.transport.faceTemperature(w.faceID(1))=tr.referenceTemperature;reject('reference-T-as-wall-T',bad);
assert(~exist(badFolder,'dir')&&strcmp(sourceSHA,fileSHA(evidenceFile))&&strcmp(legacySHA,fileSHA(legacyEvidence)));
uiAfter=cellfun(@(p)fileSHA(fullfile(root,p)),uiPaths,'UniformOutput',false);
assert(isequal(uiBefore,uiAfter),'AxisFV:SutherlandViewUI','共享UI源码不得因专用研究页而改变。');
report=struct('status','passed-actual-Sutherland-TMR-view-and-constant-JSON-regression', ...
    'sourceEvidenceSHA256',sourceSHA,'legacyEvidenceSHA256',legacySHA,'legacyJSONSHA256',fileSHA(legacyJSON), ...
    'originalRefusals',directed.refusals,'additionalRefusals',refusals,'constantRefusals',constant.refusals, ...
    'view',directed.view,'constantReplayView',constant.view, ...
    'constantJSONUnchangedExceptCurrentOrigin',true, ...
    'constantComparisonExcludedFields',{{'researchIdentity.source','researchIdentity.currentSourceBindings'}}, ...
    'constantHTMLByteIdentical',true,'sharedUISourceUnchanged',true, ...
    'referenceViscosityPaS',gas.mu,'referenceTemperatureK',tr.referenceTemperature, ...
    'wallMuRangePaS',[min(w.viscosity),max(w.viscosity)], ...
    'wallKRangeWmK',[min(w.thermalConductivity),max(w.thermalConductivity)], ...
    'sourceBindings',currentBindings(root),'browserVisuallyInspected',false,'webGLRuntimeTested',false, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'tmr_sutherland_view_evidence.mat'),'report','sourceSHA','legacySHA','uiBefore','uiAfter');
fid=fopen(fullfile(folder,'tmr_sutherland_view_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Sutherland view passed: reference/local transport distinct, 18 original + %d transport refusals, constant JSON/HTML replay, no browser/production claim.\n',numel(refusals));
    function reject(label,current)
        try
            exportTmrAxisymmetricStudy(mesh,current,saved.boundary,badFolder);
        catch exception
            assert(startsWith(exception.identifier,'AxisFV:TmrView'),'AxisFV:SutherlandViewReject', ...
                '输运篡改%s出现非预期异常%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFV:SutherlandViewNoRefusal','损坏输运%s未被拒绝。',label);
    end
end
function path=pick(candidates)
for j=1:numel(candidates),if exist(candidates{j},'file')==2,path=candidates{j};return;end,end
error('AxisFV:SutherlandViewEvidence','未找到需要的冻结证据：%s。',strjoin(candidates,' / '));
end
function data=withoutCurrentOrigin(data)
data.researchIdentity=rmfield(data.researchIdentity,{'source','currentSourceBindings'});
end
function out=currentBindings(root)
paths={'tests/testTmrAxisymmetricSutherlandView.m','tests/testTmrAxisymmetricSteadyView.m', ...
    'regionalFlow/axisymmetric/exportTmrAxisymmetricStudy.m','regionalFlow/axisymmetric/axisymmetricFVTransport.m', ...
    'regionalFlow/axisymmetric/createTmrAxisymmetricAirGas.m'};out=struct([]);
for j=1:numel(paths)
    item=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));
    if isempty(out),out=item;else,out(end+1)=item;end %#ok<AGROW>
end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
