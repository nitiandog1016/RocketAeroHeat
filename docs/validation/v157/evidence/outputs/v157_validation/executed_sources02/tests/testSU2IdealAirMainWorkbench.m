function report=testSU2IdealAirMainWorkbench(out,mode)
%TESTSU2IDEALAIRMAINWORKBENCH 主五步骤变比热三维启动、定量展示和原生取消。
% 模块一。严格拒绝旧低质量小例及错误边界，新解析壳通过原完整域门槛。
% 模块二。正常任务仅4更新，review仅复用已有数值，取消任务只启动一次。
% 模块三。核对原面、方法、温度身份、固定色标、周向导出及状态失效。
if nargin<2,mode='run';end
assert(~isfolder(out));mkdir(out);
meshFolder='outputs/v157_validation/fixture';
if isfolder('outputs/v157_validation/fixture_workbench'),meshFolder='outputs/v157_validation/fixture_workbench';end
data=load(fullfile(meshFolder,'mesh_input.mat'));surface=struct('V',data.V,'F',data.F);
base=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=base.cfg;
cfg.output.verbose=false;cfg.flight.Mach=5;cfg.flight.alpha=20;cfg.flight.beta=0;cfg.flight.Twall=300;
atm=standardAtmosphere(cfg.flight.altitude);atm.temperature=226.5;atm.pressure=1171.8;
atm.density=atm.pressure/(atm.R*atm.temperature);atm.speedOfSound=sqrt(atm.gamma*atm.R*atm.temperature);
cfg.flight.atmosphere=atm;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',4,'outputEvery',2,'namedCheckpoints',true, ...
    'maximumSeconds',60,'checkpointLimitMiB',32,'pythonExecutable','D:/anacondapython/python.exe', ...
    'initialThicknessM',.001,'turbulenceIntensity',.01,'turbulentViscosityRatio',10,'turbulentPrandtl',.9, ...
    'turbulenceProvenance','v157 explicit fully turbulent interface assumption');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(surface);
fig.UserData.actions.readHybridMesh('outputs/v145_validation/near_wall_fixture');
fig.UserData.actions.selectStep(3);ctrl=fig.UserData.controls;
ctrl.thermodynamics.Value='idealAir';ctrl.thermodynamics.ValueChangedFcn(ctrl.thermodynamics,[]);
fig.UserData.actions.confirmBoundary();
assert(~fig.UserData.workflowState.boundaryConfirmed && contains(fig.UserData.status.Text,'完整域质量未通过'));
fig.UserData.actions.readHybridMesh(meshFolder);
fig.UserData.actions.selectStep(3);ctrl=fig.UserData.controls;
ctrl.thermodynamics.Value='idealAir';ctrl.thermodynamics.ValueChangedFcn(ctrl.thermodynamics,[]);
assert(strcmp(ctrl.flowModel.Value,'SST-2003m') && strcmp(ctrl.initialMode.Value,'sstWallDistance'));
ctrl.wall.Value=2100;ctrl.wall.ValueChangedFcn(ctrl.wall,[]);fig.UserData.actions.confirmBoundary();
assert(~fig.UserData.workflowState.boundaryConfirmed);
ctrl.wall.Value=300;ctrl.wall.ValueChangedFcn(ctrl.wall,[]);
ctrl.thermalMode.Value='adiabatic';ctrl.thermalMode.ValueChangedFcn(ctrl.thermalMode,[]);
fig.UserData.actions.confirmBoundary();assert(~fig.UserData.workflowState.boundaryConfirmed);
ctrl.thermalMode.Value='isothermal';ctrl.thermalMode.ValueChangedFcn(ctrl.thermalMode,[]);
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fs=fig.UserData.workflowState.freestream;
air=idealAirThermodynamicState('PT',1171.8,226.5);assert(abs(fs.speed-5*air.sound)<2e-10);
if strcmp(mode,'run')
    fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
    folder=fig.UserData.engineJob.folder;write(fullfile(out,'engine_folder.json'),struct('folder',folder));
    finish(fig,180);assert(strcmp(fig.UserData.workflowState.status,'research-ready'),fig.UserData.status.Text);
else
    source=jsondecode(fileread('outputs/v157_validation/main02/engine_folder.json'));folder=source.folder;
    fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
end
r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;
assert(a.checkpointIdentityConfirmed,a.reason);assert(a.checkpointIteration==3 && a.cumulativeCheckpointIteration==3);
assert(r.thermalStateAudit.interfacePassed && r.thermalStateAudit.successfulUpdates==4 && ~r.productionQualified);
assert(isequal(r.V,surface.V) && isequal(r.F,surface.F) && numel(r.pressureRaw)==8);
ui=fig.UserData.researchUI;ui.setSectionX(.25);waitGL(ui);ui.inspectFace(2,[]);
assert(contains(strjoin(ui.probe.Value,newline),'OTH') && contains(ui.auditText.Text,'本段完成 4'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.low.Value=-1e5;ui.high.Value=1e5;ui.apply.ButtonPushedFcn(ui.apply,[]);waitGL(ui);
assert(isequal(ui.geometry.html.UserData.payload.faceField.range,[-1e5 1e5]));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));section=ui.getSection();
assert(size(section.segments,1)==4);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
save(fullfile(out,'result.mat'),'r','fs','-v7');
% 模块四。正常数值终态已存盘，后续读取或取消错误均不得补跑正常任务。
fig.UserData.actions.selectStep(4);fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));waitGL(fig.UserData.researchUI);
fig.UserData.actions.selectStep(3);ctrl=fig.UserData.controls;ctrl.wall.Value=301;ctrl.wall.ValueChangedFcn(ctrl.wall,[]);
assert(isempty(fig.UserData.workflowState.researchResult) && ~fig.UserData.workflowState.boundaryConfirmed);
rejected=false;try,fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));catch,rejected=true;end
assert(rejected);ctrl.wall.Value=300;ctrl.wall.ValueChangedFcn(ctrl.wall,[]);fig.UserData.actions.confirmBoundary();
fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
cancel=struct('attempted',false);
if strcmp(mode,'run') || strcmp(mode,'cancel')
    fig.UserData.actions.selectStep(3);ctrl=fig.UserData.controls;
    ctrl.iterations.Value=1000;ctrl.iterations.ValueChangedFcn(ctrl.iterations,[]);
    ctrl.outputEvery.Value=20;ctrl.outputEvery.ValueChangedFcn(ctrl.outputEvery,[]);fig.UserData.actions.confirmBoundary();
    assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
    fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
    cancelFolder=fig.UserData.engineJob.folder;write(fullfile(out,'cancel_folder.json'),struct('folder',cancelFolder));
    t=tic;requested=false;
    while fig.UserData.workflowState.busy && toc(t)<160
        drawnow;
        if isfile(fullfile(cancelFolder,'restart_flow_000000.csv')) && isfile(fullfile(cancelFolder,'surface_flow_000000.vtk'))
            receipt=jsondecode(fileread(fullfile(cancelFolder,'engine_process.json')));
            if strcmp(receipt.status,'running')
                fig.UserData.actions.cancel();requested=true;break;
            end
        end
        pause(.01);
    end
    finish(fig,90);
    native=jsondecode(fileread(fullfile(cancelFolder,'engine_process.json')));
    pipeline=jsondecode(fileread(fullfile(cancelFolder,'workflow','pipeline_process.json')));
    cancel=struct('attempted',true,'requestedWhileRunning',requested,'folder',cancelFolder,'native',native,'pipeline',pipeline);
    write(fullfile(out,'cancel.json'),cancel);
    assert(strcmp(native.status,'terminal') && strcmp(pipeline.status,'terminal'));
    assert(requested && native.cancelled && pipeline.cancelled && pipeline.exitCode==125);
    assert(~isempty(native.checkpoints.complete) && isempty(fig.UserData.workflowState.researchResult));
end
report=struct('passed',true,'mode',mode,'engineFolder',folder,'originalFaces',8,'sectionSegments',4, ...
    'actualFreestreamSpeed',fs.speed,'newModelMainFiveSteps',true,'qualityGateRejected',true, ...
    'badBoundaryRejected',true,'fixedColorScale',true,'unknownSolvedTw',true,'sameConditionReopened',true, ...
    'upstreamInvalidation',true,'cancellation',cancel,'nativeMouseObserved',false,'physicsQualified',false,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function finish(fig,seconds)
t=tic;while fig.UserData.workflowState.busy && toc(t)<seconds,drawnow;pause(.1);end
if fig.UserData.workflowState.busy,fig.UserData.actions.cancel();end
assert(~fig.UserData.workflowState.busy,'任务仍活动，已请求取消，等待实际终态。');
end
function waitGL(ui)
t=tic;while toc(t)<30
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && ~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==8 && s.rendered.webglError==0 && ~isempty(ui.getSection()));
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
