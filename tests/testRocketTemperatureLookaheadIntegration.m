function report=testRocketTemperatureLookaheadIntegration(folder)
%TESTROCKETTEMPERATURELOOKAHEADINTEGRATION 回放旧完整默认输出并检查无原态面的受限拒绝。
% 模块一：两种CV排列下与v89原函数逐值比较全部非计时对象，缺省不启用前瞻。
% 模块二：静止均匀气体的内部采样均非region1/4，原型必须明确拒绝而不漏选单元。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'),fullfile(root,'tests','fixtures'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
s=load(rocketTemperatureResearchFixture('temperature_active_unit.mat','outputs/v89_validation/active_unit01/temperature_active_unit.mat'),'cases');cases=cell(2,1);
for j=1:2
    c=s.cases{j};rhs=c.source.base.rhs;
    [old,oldReport,oldEvidence]=rocketTemperatureActiveV89Reference(c.packed,c.factor,rhs,c.source.metric,c.rows,struct('maximumIterations',900));
    [value,detail,evidence]=solveRocketTemperatureActiveDirection(c.packed,c.factor,rhs,c.source.metric,c.rows,struct('maximumIterations',900));
    assert(isequaln(old,value)&&isequaln(scrub(oldReport),scrub(detail))&&isequaln(scrub(oldEvidence),scrub(evidence)));
    assert(isequaln(value,c.value),'TemperatureLookahead:Legacy','旧保存小域结果未逐值回放。');
    cases{j}=struct('value',value,'report',detail,'evidence',evidence,'legacyValue',old,'legacyReport',oldReport,'legacyEvidence',oldEvidence);
end
a=load(rocketTemperatureResearchFixture('temperature_lookahead_unit.mat','outputs/v90_validation/unit01/temperature_lookahead_unit.mat'),'cases');small=a.cases{1};
n=size(small.Q,1);rho=.1*ones(n,1);cv=small.gas.R/(small.gas.gamma-1);
Q=[rho zeros(n,3) rho*cv*100];cut=small.cut;gas=small.gas;scale=small.scale;
ids=find(cut.openFaces.area>0&cut.openFaces.cellMinus>0&cut.openFaces.cellPlus>0);
[~,audit]=cutCellExactEulerFlux(Q(cut.openFaces.cellMinus(ids),:),Q(cut.openFaces.cellPlus(ids),:),cut.openFaces.normal(ids,:),gas);audit.faceID=ids;
assert(~any(audit.sampleRegion==1|audit.sampleRegion==4));value=zeros(n,5);value(1,5)=-rho(1)*cv*200/scale(5);
rejected=false;identifier='';try
    selectRocketTemperatureLookaheadRows(Q,struct('exactEulerFluxAudit',audit),cut,gas,scale,value,1/8);
catch e,identifier=e.identifier;rejected=strcmp(identifier,'TemperatureLookahead:PureSamples');end
assert(rejected,'TemperatureLookahead:Applicability','缺少原态面时必须明确拒绝，不删除所需单元。');
report=struct('status','passed-legacy-full-nontiming-replay-and-stationary-non-pure-sample-rejection', ...
    'allLegacyNontimingFieldsExactlyReplayed',true,'unsupportedPureSampleIdentityRejected',rejected, ...
    'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_lookahead_integration.mat'),'cases','report','Q','gas','scale','value','audit','identifier','-v7');disp(report);
end
function y=scrub(x)
y=x;
if isstruct(x)
    keys=intersect(fieldnames(x),{'elapsedSeconds','constructionSeconds'});if ~isempty(keys),y=rmfield(y,keys);end
    for j=1:numel(y),for key=fieldnames(y)',y(j).(key{1})=scrub(y(j).(key{1}));end,end
elseif iscell(x)
    for j=1:numel(x),y{j}=scrub(x{j});end
end
end
