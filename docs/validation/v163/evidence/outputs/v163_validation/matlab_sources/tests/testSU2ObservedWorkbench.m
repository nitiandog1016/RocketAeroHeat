function report=testSU2ObservedWorkbench(out)
%TESTSU2OBSERVEDWORKBENCH 新观察包沿主五步骤只读查看，复核旧保存态与异常拒绝。
% 模块一。复用v157的已验证几何/输入，避免增加数值任务。
% 模块二。核对实际包方法、原面导出及未知固体壁温，失败任务保持拒绝。
assert(~isfolder(out));mkdir(out);
source=jsondecode(fileread('outputs/v157_validation/main03/engine_folder.json'));
d=load(fullfile(source.folder,'input_snapshot.mat'),'frozen');
folder=fullfile(pwd,'outputs','v163_validation','small02');
fig=rocketAeroHeatWorkbench(d.frozen.cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.frozen.mesh);
fig.UserData.actions.readHybridMesh('outputs/v157_validation/fixture_workbench');
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.selectStep(4);fig.UserData.actions.readSSTStudy(folder);
r=fig.UserData.workflowState.researchResult;
assert(r.thermalStateAudit.interfacePassed && r.thermalStateAudit.successfulUpdates==4);
assert(strcmp(r.receipt.runtimeFlavor,'rah-ideal-air-observed-v1') && ~r.productionQualified);
assert(isequal(r.V,d.frozen.mesh.V) && isequal(r.F,d.frozen.mesh.F) && numel(r.pressureRaw)==8);
ui=fig.UserData.researchUI;ui.setSectionX(.25);waitGL(ui);ui.inspectFace(2,[]);
assert(contains(strjoin(ui.probe.Value,newline),'17位有效数字ASCII'));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
section=ui.getSection();assert(size(section.segments,1)==4);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
volume=readSU2HybridVolume('outputs/v157_validation/fixture_workbench',d.frozen.mesh);
rejected=false;
try
    readSU2VerifiedSSTStudy(fullfile(pwd,'outputs','v163_validation','real02'),fig.UserData.workflowState,volume);
catch err
    rejected=contains(err.message,'原引擎未成功结束');
end
assert(rejected);
% 模块三。旧phase-v2原始审核继续通过，历史异常首帧保持其拒绝身份。
oldReceipt=jsondecode(fileread(fullfile(source.folder,'engine_process.json')));
[~,~,old]=readSU2IdealAirStateAudit(source.folder,oldReceipt);
assert(old.interfacePassed && old.successfulUpdates==4);
assert(~isfile('outputs/v161_validation/control/ideal_air_review.json'));
save(fullfile(out,'result.mat'),'r','-v7');
report=struct('passed',true,'originalFaces',8,'sectionSegments',4,'newRuntimeMethod',true, ...
    'unknownSolvedTw',true,'abnormalTaskRejected',rejected,'legacyPhaseV2Read',true, ...
    'nativeRuns',0,'nativeMouseObserved',false,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);fileGuard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
end

function waitGL(ui)
t=tic;
while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && ~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==8 && s.rendered.webglError==0);
end

function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
