function report=testTmrAxisymmetricDurableRestart(folder,finalEvidenceFile,earlierCheckpointFile)
%TESTTMRAXISYMMETRICDURABLERESTART 在真实链副本上检查零推进重启、孤儿恢复和拒绝。
% 模块一：只接受真正收敛的完整链及它的真实早期前缀，不制造终态或放宽三步门槛。
% 模块二：原字节复制latest/全部step至私有临时目录；原both、父源和process记录不改。
% 模块三：稳态重启不调用BE；早期latest恢复真实孤儿端点至稳态，逐值还原全链结果。
% 模块四：坏源绑定、暖源SHA、step SHA/缺失、初态断接及坏孤儿必须明确拒绝。
% 外层runCount/诊断调用/真实计时/运行PID可不同；它们不能冒充新增接受步或solver成本。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'tests'),fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v67_validation','durable_restart_test');end
if nargin<2
    finalEvidenceFile=fullfile(root,'outputs','v67_validation','both','tmr_sensitivity_evidence.mat');
    if ~exist(finalEvidenceFile,'file')
        finalEvidenceFile=fullfile(root,'validation','v67','both','tmr_sensitivity_evidence.mat');
    end
end
if nargin<3
    earlierCheckpointFile=fullfile(root,'outputs','v67_validation','review','durable_checkpoint_after_segment002.mat');
    if ~exist(earlierCheckpointFile,'file')
        earlierCheckpointFile=fullfile(root,'validation','v67','review','durable_checkpoint_after_segment002.mat');
    end
end
assert(exist(finalEvidenceFile,'file')==2&&exist(earlierCheckpointFile,'file')==2, ...
    'AxisDurableTest:Evidence','需要真实最终稳态和真实早期原字节latest，不能伪造检查点。');
sourceSHA=fileSHA(finalEvidenceFile);earlierSHA=fileSHA(earlierCheckpointFile);
final=load(finalEvidenceFile);earlier=load(earlierCheckpointFile);
assert(all(isfield(final,{'report','mesh','gas','boundary','Q','initialInfo','implicit','durable'}))&& ...
    all(isfield(earlier,{'mesh','gas','boundary','Q','initialInfo','implicit','durable'})), ...
    'AxisDurableTest:Evidence','输入不是完整持久化检查点。');
requireActualSteady(final);
assert(~earlier.implicit.steadyConverged&&earlier.implicit.acceptedSteps>=1&& ...
    earlier.implicit.acceptedSteps<final.implicit.acceptedSteps&& ...
    isequaln(earlier.mesh,final.mesh)&&isequaln(earlier.gas,final.gas)&& ...
    isequaln(earlier.boundary,final.boundary)&&isequal(earlier.Q,final.Q)&& ...
    isequaln(earlier.initialInfo,final.initialInfo)&& ...
    isequaln(earlier.durable.sourceBindings,final.durable.sourceBindings)&& ...
    isequaln(earlier.durable.warmBinding,final.durable.warmBinding), ...
    'AxisDurableTest:Prefix','早期源须是同一真实链尚未收敛的前缀。');
n=final.implicit.acceptedSteps;n0=earlier.implicit.acceptedSteps;
assert(numel(final.durable.steps)==n&&numel(earlier.durable.steps)==n0&& ...
    isequaln(earlier.durable.steps,final.durable.steps(1:n0)), ...
    'AxisDurableTest:Prefix','早期step清单不是最终不可变链的逐值前缀。');
for key=fieldnames(earlier.implicit.history)'
    name=key{1};a=earlier.implicit.history.(name);b=final.implicit.history.(name);
    if strcmp(name,'stepDefect'),b=b(1:n0+1,:);else,b=b(1:n0+1);end
    assert(isequaln(a,b),'AxisDurableTest:Prefix','早期历史字段%s不是原链前缀。',name);
end
bindings=final.durable.sourceBindings;
for k=1:numel(bindings)
    assert(strcmp(fileSHA(fullfile(root,bindings(k).path)),bindings(k).sha256), ...
        'AxisDurableTest:Source','真实链的实际源码已变，不能换实现后沿用旧成功声明。');
end
warmSHA=fileSHA(final.durable.warmBinding.path);
assert(strcmp(warmSHA,final.durable.warmBinding.sha256));
sourceFolder=fileparts(finalEvidenceFile);stepHashes=cell(n,1);
for k=1:n
    relative=final.durable.steps(k).relativePath;checkStepPath(relative,k);
    stepHashes{k}=fileSHA(fullfile(sourceFolder,strrep(relative,'/',filesep)));
    assert(strcmp(stepHashes{k},final.durable.steps(k).sha256));
end

temporary=[tempname '_RocketAeroHeat_DurableRestartTest'];mkdir(temporary);
cleanup=onCleanup(@()removeTemporary(temporary)); %#ok<NASGU>
steadyFolder=fullfile(temporary,'steady_restart');orphanFolder=fullfile(temporary,'orphan_restart');
copyFixture(finalEvidenceFile,sourceFolder,final.durable.steps,steadyFolder);
copyFixture(earlierCheckpointFile,sourceFolder,final.durable.steps,orphanFolder);
refusals=struct('label',{},'identifier',{});

%% 一、已经收敛的真实latest重启：空间重放可以增加，接受步及solver成本不能增加
steadyReport=studyTmrAxisymmetricDurable(steadyFolder,1);
steady=load(fullfile(steadyFolder,'tmr_sensitivity_evidence.mat'));
assert(steadyReport.steadyConverged&&isequaln(steady.implicit,final.implicit)&& ...
    isequaln(steady.durable.steps,final.durable.steps), ...
    'AxisDurableTest:SteadyRestart','稳态重启修改Q、原history、接受步、BE账或solver计数。');
assertNoNewStepFiles(steadyFolder,n,final.durable.steps);
steadyLedger=replayAndLedger(steady);
steadyDiagnosticAdded=steady.durable.diagnosticResidualEvaluations-final.durable.diagnosticResidualEvaluations;
assert(steadyDiagnosticAdded>=1,'AxisDurableTest:Diagnostic','稳态重启的实际空间重放没有记诊断调用。');

%% 二、真实早期latest加最终原step：恢复孤儿至原三步稳态，不能执行新的BE
orphanReport=studyTmrAxisymmetricDurable(orphanFolder,1);
recovered=load(fullfile(orphanFolder,'tmr_sensitivity_evidence.mat'));
assert(orphanReport.steadyConverged&&isequaln(recovered.implicit,final.implicit)&& ...
    isequaln(recovered.durable.steps,final.durable.steps)&& ...
    recovered.durable.recoveredOrphanSteps-earlier.durable.recoveredOrphanSteps==n-n0, ...
    'AxisDurableTest:OrphanRecovery','真实孤儿恢复没有逐值还原终态/全history/solver成本，或混入新BE步。');
assertNoNewStepFiles(orphanFolder,n,final.durable.steps);
orphanLedger=replayAndLedger(recovered);
originalRecoveredDiagnostics=0;
for k=n0+1:n
    s=load(fullfile(sourceFolder,strrep(final.durable.steps(k).relativePath,'/',filesep)),'stepContext');
    originalRecoveredDiagnostics=originalRecoveredDiagnostics+s.stepContext.diagnosticResidualEvaluations;
end
orphanDiagnosticAdded=recovered.durable.diagnosticResidualEvaluations-earlier.durable.diagnosticResidualEvaluations;
% 每个孤儿的初末验证各一次，旧步已发生的诊断按原context补回；新重启也有一次重放。
expectedOrphanDiagnostics=1+2*(n-n0)+originalRecoveredDiagnostics;
assert(orphanDiagnosticAdded==expectedOrphanDiagnostics, ...
    'AxisDurableTest:Diagnostic','孤儿恢复的新增验证或补记原步诊断调用不完整。');

%% 三、同一私有夹具上的损坏输入拒绝；每次重置latest都来自真实原字节
resetLatest(steadyFolder,finalEvidenceFile);p=load(latest(steadyFolder));
p.durable.sourceBindings(1).sha256=changedHash(p.durable.sourceBindings(1).sha256);savePayload(latest(steadyFolder),p);
reject('changed-current-source-binding',@()studyTmrAxisymmetricDurable(steadyFolder,1),{'AxisDurable:Checkpoint'});

resetLatest(steadyFolder,finalEvidenceFile);p=load(latest(steadyFolder));
p.durable.warmBinding.sha256=changedHash(p.durable.warmBinding.sha256);savePayload(latest(steadyFolder),p);
reject('changed-warm-source-SHA',@()studyTmrAxisymmetricDurable(steadyFolder,1),{'AxisDurable:Checkpoint'});

resetLatest(steadyFolder,finalEvidenceFile);first=fullfile(steadyFolder,strrep(final.durable.steps(1).relativePath,'/',filesep));
fid=fopen(first,'ab');assert(fid>=0);fwrite(fid,uint8(0),'uint8');fclose(fid);
reject('changed-immutable-step-bytes',@()studyTmrAxisymmetricDurable(steadyFolder,1),{'AxisDurable:StepSHA'});
copyfile(fullfile(sourceFolder,strrep(final.durable.steps(1).relativePath,'/',filesep)),first,'f');

resetLatest(steadyFolder,finalEvidenceFile);delete(first);
reject('missing-required-step-file',@()studyTmrAxisymmetricDurable(steadyFolder,1),{'AxisDurable:Source'});
copyfile(fullfile(sourceFolder,strrep(final.durable.steps(1).relativePath,'/',filesep)),first,'f');

resetLatest(steadyFolder,finalEvidenceFile);p=load(latest(steadyFolder));p.Q(1,1)=p.Q(1,1)+1e-6;
savePayload(latest(steadyFolder),p);
reject('initial-Q-vs-initialQ-disconnected',@()studyTmrAxisymmetricDurable(steadyFolder,1),{'AxisDurable:Checkpoint'});

resetLatest(steadyFolder,finalEvidenceFile);p=load(latest(steadyFolder));p=rmfield(p,'durable');savePayload(latest(steadyFolder),p);
reject('existing-target-without-durable-identity',@()studyTmrAxisymmetricDurable(steadyFolder,1),{'AxisDurable:Checkpoint'});

%% 四、未被latest清单收录的孤儿也必须重新验方程，不能凭文件名接入
resetLatest(orphanFolder,earlierCheckpointFile);orphanPath=fullfile(orphanFolder,strrep(final.durable.steps(n0+1).relativePath,'/',filesep));
originalOrphan=fullfile(sourceFolder,strrep(final.durable.steps(n0+1).relativePath,'/',filesep));
s=load(orphanPath);s.segment.initialQ(1,1)=s.segment.initialQ(1,1)+1e-6;savePayload(orphanPath,s);
reject('orphan-initial-Q-disconnected',@()studyTmrAxisymmetricDurable(orphanFolder,1),{'AxisDurable:Orphan'});
copyfile(originalOrphan,orphanPath,'f');

resetLatest(orphanFolder,earlierCheckpointFile);s=load(orphanPath);
s.segment.Q(1,4)=s.segment.Q(1,4)+1;savePayload(orphanPath,s);
reject('orphan-stale-finalInfo-vs-Q',@()studyTmrAxisymmetricDurable(orphanFolder,1),{'AxisDurable:Orphan'});
copyfile(originalOrphan,orphanPath,'f');

resetLatest(orphanFolder,earlierCheckpointFile);s=load(orphanPath);s.stepContext.nextPseudoCFL=2*s.stepContext.nextPseudoCFL;
savePayload(orphanPath,s);
reject('orphan-next-CFL-inconsistent',@()studyTmrAxisymmetricDurable(orphanFolder,1),{'AxisDurable:Orphan'});
copyfile(originalOrphan,orphanPath,'f');

%% 五、原证据和所有step完整字节保持，测试不能给生产面积或实验资格
assert(strcmp(sourceSHA,fileSHA(finalEvidenceFile))&&strcmp(earlierSHA,fileSHA(earlierCheckpointFile))&& ...
    strcmp(warmSHA,fileSHA(final.durable.warmBinding.path)), ...
    'AxisDurableTest:SourceMutation','重启/拒绝测试修改了真实源文件。');
for k=1:n
    assert(strcmp(stepHashes{k},fileSHA(fullfile(sourceFolder,strrep(final.durable.steps(k).relativePath,'/',filesep)))));
end
report=struct('status','passed-real-steady-restart-and-real-orphan-chain-recovery', ...
    'sourceEvidenceSHA256',sourceSHA,'earlierCheckpointSHA256',earlierSHA, ...
    'sourceAcceptedSteps',n,'earlierAcceptedSteps',n0,'recoveredActualEndpoints',n-n0, ...
    'steadyRestartFullSolverResultUnchanged',true,'orphanRecoveredFullSolverResultIdentical',true, ...
    'noNewBEStepOrSolverCost',true,'steadyRestartDiagnosticCallsAdded',steadyDiagnosticAdded, ...
    'orphanDiagnosticCallsAdded',orphanDiagnosticAdded,'originalRecoveredDiagnosticCalls',originalRecoveredDiagnostics, ...
    'diagnosticRunIdentityAndElapsedNotForcedEqual',true,'steadyLedger',steadyLedger,'orphanLedger',orphanLedger, ...
    'refusals',refusals,'originalSourcesAndStepBytesUnchanged',true,'originalProcessRecordsEdited',false, ...
    'sourceBindings',sourceBindings(root),'physicalSolution',false,'experimentalValidation',false, ...
    'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'durable_restart_test_evidence.mat'),'report','sourceSHA','earlierSHA','stepHashes');
fid=fopen(fullfile(folder,'durable_restart_test_summary.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear closer;
fprintf('Durable restart passed: true %d-step prefix -> %d accepted endpoints, no new BE, %d explicit refusals; original evidence unchanged.\n', ...
    n0,n,numel(refusals));
    function reject(label,action,expected)
        try,action();catch exception
            assert(ismember(exception.identifier,expected),'AxisDurableTest:UnexpectedRefusal', ...
                '%s期望明确持久化拒绝，实际%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisDurableTest:MissingRefusal','损坏夹具%s未拒绝，测试不允许继续。',label);
    end
end

function requireActualSteady(item)
% 模块用途：不改输入flag/历史或阈值，必须已有原标准下的实际三步收敛链。
r=item.implicit;
assert(r.steadyConverged&&item.report.steadyConverged&&strcmp(r.stopReason,'steady-residual-converged')&& ...
    r.acceptedSteps>=3&&r.options.steadyResidualTolerance==2e-6&&r.options.steadyRequiredSteps==3&& ...
    r.options.minimumSteadySteps==3&&r.options.conservationTolerance==5e-10&&r.residualNorm<=2e-6&& ...
    all(r.history.residual(end-2:end)<=2e-6)&&r.globalBalanceRelative<=5e-10&& ...
    numel(r.history.residual)==r.acceptedSteps+1&&~r.physicalSolution&&~r.eligibleForSurfaceCoverage&& ...
    strcmp(item.durable.kind,'one-accepted-BE-step-checkpoint-chain-v1')&&item.durable.allAcceptedEndpointsSaved, ...
    'AxisDurableTest:NotSteady','必须等待真实链在原三步标准下收敛，不能用工作预算终态测试。');
end

function result=replayAndLedger(item)
% 模块用途：独立重放完整空间状态，并由真实初末Q及累计有单位B/S核对全五量账。
r=item.implicit;mesh=item.mesh;V=mesh.volumePerRadian;
[rate,info]=axisymmetricRANSResidual(r.Q,mesh,item.gas,item.boundary);
assert(isequaln(info,r.finalInfo)&&max(abs(rate)./r.residualScale,[],'all')==r.residualNorm, ...
    'AxisDurableTest:Replay','副本恢复后的原完整算子不能逐值重放。');
first=sum(r.initialQ.*V,1);last=sum(r.Q.*V,1);
assert(isequal(first,r.initialConservedIntegral)&&isequal(last,r.finalConservedIntegral));
defect=last-first+r.accumulatedBoundaryIntegral-r.accumulatedSourceIntegral;
relative=abs(defect)./(sum(V)*r.residualScale);
assert(all(relative<=5e-10)&&max(relative)==r.globalBalanceRelative);
result=struct('rawFiveQuantityDefect',defect,'relativeByFixedQuantityScale',relative, ...
    'fullSpatialInfoReplayed',true,'diagnosticResidualEvaluationsInThisTestHelper',1, ...
    'perRadian',true,'unitSeparated',true,'physicalTimeAdvanced',0);
end

function copyFixture(latestSource,sourceFolder,steps,target)
% 模块用途：有效夹具只复制真实原字节，保存原SHA；不从数值重建有效latest。
mkdir(target);mkdir(fullfile(target,'steps'));copyfile(latestSource,latest(target));
assert(strcmp(fileSHA(latestSource),fileSHA(latest(target))));
for k=1:numel(steps)
    rel=steps(k).relativePath;checkStepPath(rel,k);
    from=fullfile(sourceFolder,strrep(rel,'/',filesep));to=fullfile(target,strrep(rel,'/',filesep));
    copyfile(from,to);assert(strcmp(fileSHA(to),steps(k).sha256));
end
end

function assertNoNewStepFiles(folder,count,steps)
% 模块用途：原终态后不能出现新端点，恢复也不能重写任何既有不可变step。
files=dir(fullfile(folder,'steps','step_*.mat'));assert(numel(files)==count,'AxisDurableTest:NewStep','测试执行了不允许的新接受步。');
for k=1:count
    assert(strcmp(fileSHA(fullfile(folder,strrep(steps(k).relativePath,'/',filesep))),steps(k).sha256));
end
end

function checkStepPath(relative,number)
% 模块用途：只允许本链固定文件名，防止临时夹具复制/删除逃出steps目录。
assert(strcmp(relative,sprintf('steps/step_%06d.mat',number)), ...
    'AxisDurableTest:StepPath','step相对路径或顺序不符合真实链。');
end

function path=latest(folder)
path=fullfile(folder,'tmr_sensitivity_evidence.mat');
end
function resetLatest(folder,source)
copyfile(source,latest(folder),'f');assert(strcmp(fileSHA(source),fileSHA(latest(folder))));
end
function savePayload(path,payload)
% 模块用途：只在本测试临时副本中制作明确损坏负控，不覆盖真实输入。
save(path,'-struct','payload');
end
function value=changedHash(value)
if value(1)=='0',value(1)='1';else,value(1)='0';end
end
function records=sourceBindings(root)
paths={'tests/testTmrAxisymmetricDurableRestart.m','tests/studyTmrAxisymmetricDurable.m'};
records=struct('path',{},'sha256',{});
for k=1:numel(paths),records(k)=struct('path',paths{k},'sha256',fileSHA(fullfile(root,paths{k})));end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0,'AxisDurableTest:Source','缺少真实源文件。');cleanup=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');clear cleanup;md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
function removeTemporary(path)
% 模块用途：只清理已经规范化确认属于系统临时目录的本测试私有副本。
absolute=char(java.io.File(path).getCanonicalPath());temporary=char(java.io.File(tempdir).getCanonicalPath());
if startsWith(lower(absolute),[lower(temporary) filesep])&&endsWith(absolute,'_RocketAeroHeat_DurableRestartTest')&&exist(absolute,'dir')
    rmdir(absolute,'s');
end
end
