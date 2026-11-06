function report=testSU2ConservativeStepWorkbench(out)
%TESTSU2CONSERVATIVESTEPWORKBENCH 实际整箭限步包的主五步骤结果验收。
% 模块一。复用原19752面、混合体网格和工况，只读取本次正常终态。
% 模块二。核对原面云图、截面、数值方法、未知固体壁温和导出。
% 模块三。错误壁温和旧异常任务拒绝，旧phase-v2热审核继续兼容。
assert(~isfolder(out));mkdir(out);
folder=fullfile(pwd,'outputs','v165_validation','real01');
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
d=load('outputs/v118_validation/actual01/surface.mat','s');
b=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=b.cfg;
cfg.output.verbose=false;cfg.flight.Mach=5;cfg.flight.alpha=20;cfg.flight.beta=0;cfg.flight.Twall=300;
atm=standardAtmosphere(cfg.flight.altitude);atm.temperature=manifest.flight.temperature;atm.pressure=manifest.flight.pressure;
atm.density=atm.pressure/(atm.R*atm.temperature);atm.speedOfSound=sqrt(atm.gamma*atm.R*atm.temperature);
cfg.flight.atmosphere=atm;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',50,'outputEvery',50,'namedCheckpoints',true, ...
    'maximumSeconds',300,'checkpointLimitMiB',1024,'pythonExecutable','D:/anacondapython/python.exe', ...
    'thermodynamicModel','idealAir','flowModel','SST-2003m','initialMode','sstWallDistance', ...
    'initialThicknessM',.001,'turbulenceIntensity',.01,'turbulentViscosityRatio',10,'turbulentPrandtl',.9, ...
    'turbulenceProvenance','v165 saved complete domain; fully turbulent assumption; transition unknown');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v145_validation/mesh');
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.selectStep(4);fig.UserData.actions.readSSTStudy(folder);
r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;
assert(a.checkpointIdentityConfirmed,a.reason);
assert(a.checkpointIteration==49 && a.cumulativeCheckpointIteration==199 && a.completedUpdatesBeforeRun==150);
assert(r.thermalStateAudit.interfacePassed && r.thermalStateAudit.successfulUpdates==50 && ~r.productionQualified);
assert(strcmp(r.receipt.runtimeFlavor,'rah-ideal-air-step-v1'));
assert(isequal(r.V,d.s.V) && isequal(r.F,d.s.F) && numel(r.pressureRaw)==19752);
ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);ui.inspectFace(12000,[]);
probe=strjoin(ui.probe.Value,newline);
assert(contains(probe,'共同守恒方向有限半步') && contains(probe,'限步不证明全局守恒') && contains(probe,'按当前k重算') && contains(probe,'17位有效数字'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.low.Value=-1e5;ui.high.Value=1e5;ui.apply.ButtonPushedFcn(ui.apply,[]);waitGL(ui);
assert(isequal(ui.geometry.html.UserData.payload.faceField.range,[-1e5 1e5]));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
section=ui.getSection();assert(size(section.segments,1)==92);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
save(fullfile(out,'result.mat'),'r','-v7');
volume=readSU2HybridVolume('outputs/v145_validation/mesh',d.s);
bad=fig.UserData.workflowState;bad.cfg.flight.Twall=301;wrongWall=false;
try,readSU2VerifiedSSTStudy(folder,bad,volume);catch,wrongWall=true;end
assert(wrongWall);
oldFailure=false;
try
    readSU2VerifiedSSTStudy(fullfile(pwd,'outputs','v163_validation','real02'),fig.UserData.workflowState,volume);
catch err
    oldFailure=contains(err.message,'原引擎未成功结束');
end
assert(oldFailure);
source=jsondecode(fileread('outputs/v157_validation/main03/engine_folder.json'));
oldReceipt=jsondecode(fileread(fullfile(source.folder,'engine_process.json')));
[~,~,old]=readSU2IdealAirStateAudit(source.folder,oldReceipt);
assert(old.interfacePassed && old.successfulUpdates==4);
report=struct('passed',true,'originalFaces',19752,'sectionSegments',92,'checkpoint',49, ...
    'updates',50,'priorUpdates',150,'cumulativeUpdates',200,'newNumericalMethod',true, ...
    'fixedColorScale',true,'unknownSolvedTw',true,'wrongWallRejected',wrongWall, ...
    'oldAbnormalTaskRejected',oldFailure,'legacyPhaseV2Read',true, ...
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
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==19752 && s.rendered.webglError==0);
end

function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
