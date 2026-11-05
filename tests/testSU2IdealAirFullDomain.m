function report=testSU2IdealAirFullDomain(out)
%TESTSU2IDEALAIRFULLDOMAIN 真实完整箭体变比热结果的主工作台只读验收。
% 模块一。复用真实原面、一般混合网格和已保存任务，明确实际物性及边界。
% 模块二。核对全部原面、周向曲线、探针、固定色标及未知固体壁温。
% 模块三。拒绝错误工况/模型，回归本次修改影响的旧定比热只读入口。
assert(~isfolder(out));mkdir(out);
folder=fullfile(pwd,'outputs','v160_validation','real01');
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
d=load('outputs/v118_validation/actual01/surface.mat','s');
volume=readSU2HybridVolume('outputs/v145_validation/mesh',d.s);
capacity=estimateSU2CheckpointStorage(volume,150,75,true,1536,'SST-2003m','idealAir');
assert(capacity.estimatedBytes==1206232576 && capacity.withinLimit);
smallBudget=estimateSU2CheckpointStorage(volume,150,75,true,1024,'SST-2003m','idealAir');assert(~smallBudget.withinLimit);
rejected=false;try,estimateSU2CheckpointStorage(volume,150,75,true,1537,'SST-2003m','idealAir');catch,rejected=true;end
assert(rejected);
b=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=b.cfg;
cfg.output.verbose=false;cfg.flight.Mach=5;cfg.flight.alpha=20;cfg.flight.beta=0;cfg.flight.Twall=300;
atm=standardAtmosphere(cfg.flight.altitude);atm.temperature=manifest.flight.temperature;atm.pressure=manifest.flight.pressure;
atm.density=atm.pressure/(atm.R*atm.temperature);atm.speedOfSound=sqrt(atm.gamma*atm.R*atm.temperature);
cfg.flight.atmosphere=atm;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',150,'outputEvery',75,'namedCheckpoints',true, ...
    'maximumSeconds',600,'checkpointLimitMiB',1536,'pythonExecutable','D:/anacondapython/python.exe', ...
    'thermodynamicModel','idealAir','flowModel','SST-2003m','initialMode','sstWallDistance', ...
    'initialThicknessM',.001,'turbulenceIntensity',.01,'turbulentViscosityRatio',10,'turbulentPrandtl',.9, ...
    'turbulenceProvenance','v160 existing fully turbulent full-domain state; transition unknown');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v145_validation/mesh');
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.selectStep(4);fig.UserData.actions.readSSTStudy(folder);
r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;
assert(a.checkpointIdentityConfirmed,a.reason);
assert(a.checkpointIteration==149 && a.cumulativeCheckpointIteration==149 && a.completedUpdatesBeforeRun==0);
assert(r.thermalStateAudit.interfacePassed && r.thermalStateAudit.successfulUpdates==150 && ~r.productionQualified);
assert(isequal(r.V,d.s.V) && isequal(r.F,d.s.F) && numel(r.pressureRaw)==19752);
ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);ui.inspectFace(12000,[]);
probe=strjoin(ui.probe.Value,newline);assert(contains(probe,'SST-2003m') && contains(probe,'OTH'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.low.Value=-1e5;ui.high.Value=1e5;ui.apply.ButtonPushedFcn(ui.apply,[]);waitGL(ui);
assert(isequal(ui.geometry.html.UserData.payload.faceField.range,[-1e5 1e5]));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
section=ui.getSection();assert(size(section.segments,1)==92);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
save(fullfile(out,'result.mat'),'r','-v7');
bad=fig.UserData.workflowState;bad.cfg.flight.Twall=301;wrongWall=false;
try,readSU2VerifiedSSTStudy(folder,bad,volume);catch,wrongWall=true;end
assert(wrongWall);
bad=fig.UserData.workflowState;bad.cfg.su2.thermodynamicModel='idealGas';wrongModel=false;
try,readSU2VerifiedSSTStudy(folder,bad,volume);catch,wrongModel=true;end
assert(wrongModel);
bad=fig.UserData.workflowState;bad.freestream.flowVelocityVector=bad.freestream.flowVelocityVector*.99;wrongSpeed=false;
try,readSU2VerifiedSSTStudy(folder,bad,volume);catch,wrongSpeed=true;end
assert(wrongSpeed);
fig.UserData.actions.selectStep(3);ctrl=fig.UserData.controls;
ctrl.thermodynamics.Value='idealGas';ctrl.thermodynamics.ValueChangedFcn(ctrl.thermodynamics,[]);
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed);
fig.UserData.actions.readSSTStudy(fullfile(pwd,'outputs','v146_validation','control'));
old=fig.UserData.workflowState.researchResult;
assert(old.convergenceAudit.checkpointIteration==99 && numel(old.pressureRaw)==19752 && ~old.productionQualified);
report=struct('passed',true,'originalFaces',19752,'sectionSegments',92,'checkpoint',149, ...
    'updates',150,'priorUpdates',0,'newModelFullDomainRead',true,'fixedColorScale',true,'unknownSolvedTw',true, ...
    'wrongWallRejected',wrongWall,'wrongModelRejected',wrongModel,'wrongSpeedRejected',wrongSpeed, ...
    'legacyConstantCpRead',true,'solverRuns',0,'nativeMouseObserved',false,'physicsQualified',false,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function waitGL(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && ~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==19752 && s.rendered.webglError==0 && ~isempty(ui.getSection()));
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
