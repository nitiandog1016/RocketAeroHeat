function report=testSU2SSTPocketWorkbench(out)
%TESTSU2SSTPOCKETWORKBENCH 将两个完整域初场方案接入同一主工作台，不重跑流场。
% 模块一：核对物理边界/原面/实际步号，初场修正身份不隐藏。
% 模块二：同一结果区的云图、探针、92段周向曲线和导出保留方案身份。
% 模块三：未知壁温保持NaN；改变壁温的只读请求拒绝，旧普通SST结果继续兼容。
assert(~isfolder(out));mkdir(out);
d=load('outputs/v118_validation/actual01/surface.mat');c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;
source=jsondecode(fileread('outputs/v126_validation/main01/engine_folder.json'));
cfg.output.verbose=false;cfg.su2=struct('volumeMode','hybridTriangles','iterations',10,'outputEvery',10, ...
    'namedCheckpoints',true,'maximumSeconds',600,'checkpointLimitMiB',950,'pythonExecutable','D:/anacondapython/python.exe', ...
    'flowModel','SST-2003m','initialMode','sstTransfer','checkpointSource',source.folder);
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v120_validation/mesh32');
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
items=struct;
for label={'control','repaired'}
    name=label{1};suffix=['pocket_' name];if strcmp(name,'repaired'),suffix=[suffix '_roundoff'];end
    folder=fullfile(pwd,'outputs','v128_validation',suffix);
    fig.UserData.actions.selectStep(4);assert(isfield(fig.UserData.controls,'readSSTStudy'));
    fig.UserData.actions.readSSTStudy(folder);r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;
    assert(a.checkpointIdentityConfirmed,a.reason);assert(a.checkpointIteration==9 && a.sourceModelLocalIteration==149);
    assert(strcmp(a.initializationStudy.variant,name) && ~isfield(a,'cumulativeCheckpointIteration'));
    assert(~r.productionQualified && numel(r.heatFluxRaw)==19752 && r.flight.iterations==10 && r.flight.outputEvery==10);
    ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);ui.inspectFace(12000,[]);
    assert(contains(strjoin(ui.probe.Value,newline),'初场对照'));
    ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
    target=fullfile(out,name);mkdir(target);ui.exportTo(fullfile(target,'fields.csv'));ui.exportSectionTo(fullfile(target,'section.csv'));
    section=ui.getSection();assert(size(section.segments,1)==92);
    ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
    assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
    save(fullfile(target,'result.mat'),'r','-v7');
    items.(name)=struct('faces',19752,'sectionSegments',92,'sourceIteration',149,'localIteration',9, ...
        'actualRunControls',r.actualRunControls,'mainWorkbenchReadOnlyAccepted',true,'physicsQualified',false);
end
volume=readSU2HybridVolume('outputs/v120_validation/mesh32',d.s);bad=fig.UserData.workflowState;bad.cfg.flight.Twall=301;rejected=false;
try,readSU2VerifiedSSTStudy(folder,bad,volume);catch,rejected=true;end
assert(rejected,'不同给定Tw不能读取为同工况。');
old=load(fullfile(source.folder,'mapped_research.mat'),'result');a=auditSU2ResultConvergence(source.folder,old.result);
assert(a.checkpointIdentityConfirmed && isfield(a,'turbulenceRestart') && ~isfield(a,'numericalStudy'));
report=struct('variants',items,'wrongWallBoundaryRejected',rejected,'oldSSTRestartCompatible',true, ...
    'solverRuns',0,'newNativeScreenshots',false,'physicsQualified',false,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function waitGL(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision&&s.rendered.faces==19752&&s.rendered.webglError==0&&~isempty(ui.getSection()));
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
