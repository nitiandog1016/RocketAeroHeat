function report=testSU2SSTStartupWorkbench(out)
%TESTSU2SSTSTARTUPWORKBENCH 在主五步骤核对实际启动方案和原面结果，不重跑流场。
% 模块一：复用已验几何及物理边界，源与本段步号、流动阶数分别显示。
% 模块二：原面p/q/相关温度、周向曲线、探针和CSV带有准确方法来源。
% 模块三：未知求解Tw不填零，错误物理边界拒绝；旧v145初场身份保持。
assert(~isfolder(out));mkdir(out);
d=load('outputs/v118_validation/actual01/surface.mat');
c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;cfg.output.verbose=false;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',100,'outputEvery',50,'namedCheckpoints',true, ...
    'maximumSeconds',480,'checkpointLimitMiB',950,'pythonExecutable','D:/anacondapython/python.exe', ...
    'flowModel','SST-2003m','initialMode','sstWallDistance','initialThicknessM',.001,'turbulenceIntensity',.01, ...
    'turbulentViscosityRatio',10,'turbulentPrandtl',.9,'turbulenceProvenance','v146 fully turbulent research assumption');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v145_validation/mesh');
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed);
items=struct;labels={'control','first_order'};
if isfile('outputs/v146_validation/restore_second_order/sst_state_audit.json'),labels{end+1}='restore_second_order';end
for j=1:numel(labels)
    name=labels{j};folder=fullfile(pwd,'outputs','v146_validation',name);
    fig.UserData.actions.selectStep(4);fig.UserData.actions.readSSTStudy(folder);
    r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;assert(a.checkpointIdentityConfirmed,a.reason);
    assert(a.checkpointIteration==99 && ~isfield(a,'cumulativeCheckpointIteration'));
    assert(strcmp(a.numericalStudy.variant,name) && ~r.productionQualified && numel(r.heatFluxRaw)==19752);
    assert(r.flight.iterations==100 && r.flight.outputEvery==50 && r.thermalStateAudit.interfacePassed);
    ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);ui.inspectFace(12000,[]);
    probe=strjoin(ui.probe.Value,newline);assert(contains(probe,'SST-2003m'));
    if strcmp(name,'first_order'),assert(contains(probe,'一阶流动启动候选'));else,assert(contains(probe,'二阶流动重构'));end
    ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
    target=fullfile(out,name);mkdir(target);ui.exportTo(fullfile(target,'fields.csv'));ui.exportSectionTo(fullfile(target,'section.csv'));
    section=ui.getSection();assert(size(section.segments,1)==92);
    ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
    assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));save(fullfile(target,'result.mat'),'r','-v7');
    items.(name)=struct('faces',19752,'sectionSegments',92,'sourceIteration',a.sourceModelLocalIteration, ...
        'localIteration',99,'flowSpatialOrder',a.numericalStudy.flowSpatialOrder,'webgl',ui.geometry.html.UserData.rendered);
end
volume=readSU2HybridVolume('outputs/v145_validation/mesh',d.s);bad=fig.UserData.workflowState;bad.cfg.flight.Twall=301;rejected=false;
try,readSU2VerifiedSSTStudy(folder,bad,volume);catch,rejected=true;end
assert(rejected);
old=jsondecode(fileread('outputs/v145_validation/main01/engine_folder.json'));saved=load(fullfile(old.folder,'mapped_research.mat'),'result');
a=auditSU2ResultConvergence(old.folder,saved.result);assert(a.checkpointIdentityConfirmed && isfield(a,'independentInitialization'));
report=struct('variants',items,'wrongWallBoundaryRejected',rejected,'oldIndependentInitializationCompatible',true, ...
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
