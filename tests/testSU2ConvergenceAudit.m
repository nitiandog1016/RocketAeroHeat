function testSU2ConvergenceAudit(out,capture)
%TESTSU2CONVERGENCEAUDIT 用真实终态及显式反例检查时刻身份、门槛和主UI。
% 模块一：复用v109超时与v108收敛平板，不重复数值运行。
% 模块二：反例只存在test_fixture目录，文件变更、稀疏历史、重复HF均不得误通过。
% 模块三：主窗口复用v107同工况快照，核对新增说明和两种导出的身份。
if nargin<2,capture=false;end
assert(~isfolder(out));mkdir(out);
s=load('outputs/v109_validation/checkpoint01/checkpoint_original_faces.mat','result');
folder='outputs/v109_validation/rocket_seed01';a=auditSU2ResultConvergence(folder,s.result);
write(fullfile(out,'timeout_audit.json'),a);
assert(strcmp(a.status,'checkpoint-and-history-audited'),a.reason);
assert(a.checkpointIdentityConfirmed && a.checkpointIteration==3750 && a.terminalIteration==3921);
assert(~a.checkpointResidualPassed && ~a.terminalResidualPassed && ~a.productionQualified);
r=s.result;r.convergenceAudit=a;lines=describeSU2Convergence(r);
assert(contains(lines{1},'超时') && contains(lines{1},'3921') && contains(lines{2},'3750'));
plateFolder='outputs/v108_validation/plate_warm01';wall=readSU2SurfaceResult(fullfile(plateFolder,'surface_flow.vtk'));
[x,id]=sort(wall.vertices(:,1));p=struct('receipt',jsondecode(fileread(fullfile(plateFolder,'engine_process.json'))), ...
    'reconstructedEngineHeatPowerRaw',trapz(x,wall.heatFluxRaw(id)));
b=auditSU2ResultConvergence(plateFolder,p);write(fullfile(out,'plate_audit.json'),b);
% 已收敛平板的HF进入平台，17个保存候选在输出精度内无法区别。
% 此时末次残差可通过，但不能用热负荷反推云图时刻；保持未知才是正确行为。
assert(strcmp(b.status,'history-audited-checkpoint-unconfirmed'),b.reason);
assert(b.terminalResidualPassed && b.terminalIteration==23477 && b.heatLoadCandidateCount>1 && ...
    ~b.checkpointIdentityConfirmed && isnan(b.checkpointIteration) && ~b.productionQualified);
missing=auditSU2ResultConvergence(out,s.result);assert(strcmp(missing.status,'evidence-unavailable') && ~missing.checkpointIdentityConfirmed);
% 测试反例为显式合成记录；不冒充另一次引擎任务，不修改真实来源。
fixture=fullfile(out,'test_fixture');mkdir(fixture);
for name={'case.cfg','history.csv','surface_flow.vtk'},copyfile(fullfile(folder,name{1}),fullfile(fixture,name{1}));end
raw=fileread(fullfile(folder,'engine_process.json'));
raw=['{"syntheticTestFixture":true,"note":"test-only mutated evidence; no solver process",' raw(2:end)];
rr=jsondecode(raw);writeText(fullfile(fixture,'engine_process.json'),raw);testResult=s.result;testResult.receipt=rr;
fid=fopen(fullfile(fixture,'surface_flow.vtk'),'a');fprintf(fid,'\n');fclose(fid);
c=auditSU2ResultConvergence(fixture,testResult);assert(strcmp(c.status,'evidence-rejected') && ~c.checkpointIdentityConfirmed);
copyfile(fullfile(folder,'surface_flow.vtk'),fullfile(fixture,'surface_flow.vtk'));
text=fileread(fullfile(fixture,'history.csv'));rows=splitlines(string(text));header=rows(1);
% 重复同一保存HF，保留每行迭代号；SHA以合成回执绑定，考察唯一性而非只考察散列。
first=split(rows(2),',');checkpoint=split(rows(3752),',');first(end-1)=checkpoint(end-1);rows(2)=join(first,',');
writeText(fullfile(fixture,'history.csv'),char(join(rows,newline)));
changed=strrep(raw,s.result.receipt.outputSHA256.history_csv,sha(fullfile(fixture,'history.csv')));
writeText(fullfile(fixture,'engine_process.json'),changed);testResult.receipt=jsondecode(changed);
d=auditSU2ResultConvergence(fixture,testResult);
assert(strcmp(d.status,'history-audited-checkpoint-unconfirmed') && d.heatLoadCandidateCount>=2 && isnan(d.checkpointIteration));
writeText(fullfile(fixture,'history.csv'),char(join([header;rows(2);rows(4)],newline)));
changed=strrep(raw,s.result.receipt.outputSHA256.history_csv,sha(fullfile(fixture,'history.csv')));
writeText(fullfile(fixture,'engine_process.json'),changed);testResult.receipt=jsondecode(changed);
e=auditSU2ResultConvergence(fixture,testResult);assert(strcmp(e.status,'evidence-rejected') && contains(e.reason,'完整连续历史'));
write(fullfile(out,'negative_cases.json'),struct('changedOutputRejected',c,'duplicateHeatLoadUnconfirmed',d,'sparseHistoryRejected',e,'missingEvidence',missing));
input=load('docs/validation/v107/replay_seed/input_snapshot.mat','frozen');
fig=rocketAeroHeatWorkbench(input.frozen.cfg);guard=onCleanup(@()close(fig));
fig.UserData.actions.loadNative();fig.UserData.actions.generateMesh();fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
fig.UserData.actions.readResearch('docs/validation/v107/replay_seed/mapped_research.mat');
ui=fig.UserData.researchUI;state=fig.UserData.workflowState;aa=state.researchResult.convergenceAudit;
assert(strcmp(aa.status,'checkpoint-and-history-audited'),aa.reason);
assert(aa.checkpointIteration==1 && ~aa.checkpointResidualPassed && isempty(fig.UserData.engineJob));
assert(contains(ui.auditText.Text,'显示检查点 1 步') && contains(ui.auditText.Text,'未获物理资格'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(.7);
start=tic;while isempty(ui.getSection()) && toc(start)<40,drawnow;pause(.1);end
assert(~isempty(ui.getSection()));ui.exportTo(fullfile(out,'current_fields.csv'));ui.exportSectionTo(fullfile(out,'current_section.csv'));
for name={'current_fields.csv.json','current_section.csv.json'}
    % JSON将NaN编码为null，回读为空数组；比较相同序列化合同，不能要求MATLAB内存形状相同。
    meta=jsondecode(fileread(fullfile(out,name{1})));assert(isequaln(meta.convergenceAudit,jsondecode(jsonencode(aa))));
end
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('name',fig.Name));fprintf('V110_AUDIT_UI_READY\n');
    start=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(start)<60,drawnow;pause(.1);end
end
write(fullfile(out,'report.json'),struct('mainUI',true,'solverRun',false,'snapshot','v107 M4 alpha20 beta5', ...
    'rendered',ui.geometry.html.UserData.rendered,'audit',aa,'auditText',ui.auditText.Text,'newQualifiedArea',0));
fprintf('Convergence identity audit and main UI passed using existing results only.\n');
end
function write(file,value)
writeText(file,[jsonencode(value) newline]);
end
function writeText(file,text)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));b=fread(fid,inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
