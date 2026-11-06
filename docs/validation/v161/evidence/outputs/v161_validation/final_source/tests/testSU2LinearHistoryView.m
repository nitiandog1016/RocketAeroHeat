function report=testSU2LinearHistoryView(out)
%TESTSU2LINEARHISTORYVIEW 核对实际线性方法展示与失败任务拒绝。
% 模块一。复用已审核v160完整场，新增方法标签随原面探针和CSV导出。
% 模块二。v161异常和热审核失败保留拒绝，旧定比热入口继续核验。
assert(~isfolder(out));mkdir(out);
folder=fullfile(pwd,'outputs','v160_validation','real01');
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
d=load('outputs/v118_validation/actual01/surface.mat','s');
b=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=b.cfg;
cfg.output.verbose=false;cfg.flight.Mach=5;cfg.flight.alpha=20;cfg.flight.beta=0;cfg.flight.Twall=300;
atm=standardAtmosphere(cfg.flight.altitude);atm.temperature=manifest.flight.temperature;atm.pressure=manifest.flight.pressure;
atm.density=atm.pressure/(atm.R*atm.temperature);atm.speedOfSound=sqrt(atm.gamma*atm.R*atm.temperature);
cfg.flight.atmosphere=atm;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',150,'outputEvery',75,'namedCheckpoints',true, ...
    'maximumSeconds',600,'checkpointLimitMiB',1536,'pythonExecutable','D:/anacondapython/python.exe', ...
    'thermodynamicModel','idealAir','flowModel','SST-2003m','initialMode','sstWallDistance', ...
    'initialThicknessM',.001,'turbulenceIntensity',.01,'turbulentViscosityRatio',10,'turbulentPrandtl',.9, ...
    'turbulenceProvenance','v161 read-only linear-method review; fully turbulent assumption');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v145_validation/mesh');
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.selectStep(4);fig.UserData.actions.readSSTStudy(folder);
r=fig.UserData.workflowState.researchResult;
assert(isequal(r.F,d.s.F) && r.convergenceAudit.checkpointIteration==149);
assert(strcmp(r.linearSolver.LINEAR_SOLVER,'BCGSTAB') && strcmp(r.linearSolver.LINEAR_SOLVER_PREC,'LU_SGS'));
ui=fig.UserData.researchUI;ui.inspectFace(12000,[]);
probe=strjoin(ui.probe.Value,newline);assert(contains(probe,'BCGSTAB') && contains(probe,'LU_SGS'));
ui.exportTo(fullfile(out,'fields.csv'));
rejected=false;
try,fig.UserData.actions.readSSTStudy(fullfile(pwd,'outputs','v161_validation','control'));catch,rejected=true;end
assert(rejected && fig.UserData.workflowState.researchResult.convergenceAudit.checkpointIteration==149);
fig.UserData.actions.selectStep(3);ctrl=fig.UserData.controls;
ctrl.thermodynamics.Value='idealGas';ctrl.thermodynamics.ValueChangedFcn(ctrl.thermodynamics,[]);
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed);
fig.UserData.actions.readSSTStudy(fullfile(pwd,'outputs','v146_validation','control'));
old=fig.UserData.workflowState.researchResult;
assert(old.convergenceAudit.checkpointIteration==99 && strcmp(old.linearSolver.LINEAR_SOLVER,'BCGSTAB'));
report=struct('passed',true,'displayedSourceVersion',160,'originalFaces',19752,'linearMethodProbe',true, ...
    'failedV161Rejected',rejected,'legacyConstantCpRead',true,'nativeMouseObserved',false,'solverRuns',0, ...
    'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
