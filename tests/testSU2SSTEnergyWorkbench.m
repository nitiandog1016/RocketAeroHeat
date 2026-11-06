function report=testSU2SSTEnergyWorkbench(out)
%TESTSU2SSTENERGYWORKBENCH SST共同能量约束的主五步骤只读验收。
% 模块一。正常整箭用原面；异常整箭明确拒绝，使用正常小例验证方法身份。
% 模块二。核对原面、固定色标、曲线、未知固体壁温和CSV来源。
% 模块三。旧同包冷重启、旧phase-v2及错误边界保持原语义。
assert(~isfolder(out));mkdir(out);
realFolder=fullfile(pwd,'outputs','v167_validation','real01');
realReceipt=jsondecode(fileread(fullfile(realFolder,'engine_process.json')));
assert(strcmp(realReceipt.status,'terminal'));
full=realReceipt.exitCode==0 && realReceipt.checkpoints.normalScheduleComplete;
source=jsondecode(fileread('outputs/v157_validation/main03/engine_folder.json'));
if full
    folder=realFolder;manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
    d=load('outputs/v118_validation/actual01/surface.mat','s');surface=d.s;
    b=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=b.cfg;
    cfg.output.verbose=false;cfg.flight.Mach=5;cfg.flight.alpha=20;cfg.flight.beta=0;cfg.flight.Twall=300;
    atm=standardAtmosphere(cfg.flight.altitude);atm.temperature=manifest.flight.temperature;atm.pressure=manifest.flight.pressure;
    atm.density=atm.pressure/(atm.R*atm.temperature);atm.speedOfSound=sqrt(atm.gamma*atm.R*atm.temperature);cfg.flight.atmosphere=atm;
    cfg.su2=struct('volumeMode','hybridTriangles','iterations',200,'outputEvery',100,'namedCheckpoints',true, ...
        'maximumSeconds',600,'checkpointLimitMiB',1536,'pythonExecutable','D:/anacondapython/python.exe', ...
        'thermodynamicModel','idealAir','flowModel','SST-2003m','initialMode','sstWallDistance', ...
        'initialThicknessM',.001,'turbulenceIntensity',.01,'turbulentViscosityRatio',10,'turbulentPrandtl',.9, ...
        'turbulenceProvenance','v167 saved candidate; fully turbulent assumption; transition unknown');
    volumeFolder='outputs/v145_validation/mesh';faceCount=19752;segments=92;station=-1.5;probeFace=12000;last=199;prior=200;
else
    folder=fullfile(pwd,'outputs','v167_validation','small01');
    d=load(fullfile(source.folder,'input_snapshot.mat'),'frozen');cfg=d.frozen.cfg;surface=d.frozen.mesh;
    volumeFolder='outputs/v157_validation/fixture_workbench';faceCount=8;segments=4;station=.25;probeFace=2;last=2;prior=4;
end
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(surface);fig.UserData.actions.readHybridMesh(volumeFolder);
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.selectStep(4);fig.UserData.actions.readSSTStudy(folder);
r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;
assert(a.checkpointIdentityConfirmed,a.reason);
assert(a.checkpointIteration==last && a.completedUpdatesBeforeRun==prior && a.cumulativeCheckpointIteration==prior+last);
assert(r.thermalStateAudit.interfacePassed && r.thermalStateAudit.successfulUpdates==last+1 && ~r.productionQualified);
assert(strcmp(r.receipt.runtimeFlavor,'rah-sst-energy-v1') && isfield(r,'sourceCheckpointIdentity') && ~isfield(r,'restartTransition'));
assert(isequal(r.V,surface.V) && isequal(r.F,surface.F) && numel(r.pressureRaw)==faceCount);
ui=fig.UserData.researchUI;ui.setSectionX(station);waitGL(ui,faceCount);ui.inspectFace(probeFace,[]);
probe=strjoin(ui.probe.Value,newline);
assert(contains(probe,'k/ω共同松弛') && contains(probe,'SST约束候选冷重启') && contains(probe,'算子改变') && contains(probe,'非无中断等价'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.low.Value=-1e5;ui.high.Value=1e5;ui.apply.ButtonPushedFcn(ui.apply,[]);waitGL(ui,faceCount);
assert(isequal(ui.geometry.html.UserData.payload.faceField.range,[-1e5 1e5]));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
section=ui.getSection();assert(size(section.segments,1)==segments);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui,faceCount);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
volume=readSU2HybridVolume(volumeFolder,surface);
bad=fig.UserData.workflowState;bad.cfg.flight.Twall=301;wrongWall=false;
try,readSU2VerifiedSSTStudy(folder,bad,volume);catch,wrongWall=true;end
assert(wrongWall);
failedFolder=fullfile(pwd,'outputs','v166_validation','real02');if ~full,failedFolder=realFolder;end
rejected=false;
try,readSU2VerifiedSSTStudy(failedFolder,fig.UserData.workflowState,volume);
catch err,rejected=contains(err.message,'原引擎未成功结束');end
assert(rejected);
if ~full
    old=readSU2VerifiedSSTStudy(fullfile(pwd,'outputs','v166_validation','small01'),fig.UserData.workflowState,volume);
    assert(old.restartTransition.unchangedOperator && old.thermalStateAudit.interfacePassed);
end
oldReceipt=jsondecode(fileread(fullfile(source.folder,'engine_process.json')));
[~,~,old]=readSU2IdealAirStateAudit(source.folder,oldReceipt);assert(old.interfacePassed);
save(fullfile(out,'result.mat'),'r','-v7');
report=struct('passed',true,'originalFaces',faceCount,'sectionSegments',segments,'fullDomainDisplayed',full, ...
    'checkpoint',last,'priorUpdates',prior,'newNumericalMethod',true,'fixedColorScale',true, ...
    'unknownSolvedTw',true,'wrongWallRejected',wrongWall,'abnormalTaskRejected',rejected, ...
    'legacySameRuntimeRead',~full,'legacyPhaseV2Read',true,'nativeMouseObserved',false,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);fileGuard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end

function waitGL(ui,faces)
t=tic;
while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && ~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==faces && s.rendered.webglError==0);
end

function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
