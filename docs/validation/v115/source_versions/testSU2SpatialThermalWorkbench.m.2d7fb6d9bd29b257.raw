function report=testSU2SpatialThermalWorkbench(out,capture)
%TESTSU2SPATIALTHERMALWORKBENCH 复用经典体网格，测试常量极限及真实空间壁温主流程。
% 模块一：官方标量与受控常量表使用同网格/工况/20步，不降低收敛门槛。
% 模块二：非均匀表由主工作台确认并实际运行，核对全部节点输入和逐面身份。
% 模块三：导出原面/截面及拾取；合成表只验接口，不是实验或物理资格。
if nargin<2,capture=false;end
assert(~isfolder(out));mkdir(out);[surface,truth]=buildSchaeferFergusonBody();condition=schaeferFergusonCondition();
meshFolder='outputs/hybrid_mesh/a75f5d87-b861-4478-b972-3095018283e6';
volume=readSU2HybridVolume(meshFolder,surface);L=truth.totalLength;
flight=struct('Mach',condition.Mach,'temperature',condition.atmosphere.temperature, ...
    'pressure',condition.atmosphere.pressure,'wallTemperature',100,'iterations',20,'outputEvery',20);
names={'official_constant','native_constant'};values=cell(2,1);receipts=cell(2,1);
for j=1:2
    if j==2,flight.thermalBoundary=struct('mode','axialTemperature','x',[0 L],'temperature',[100 100],'source','常量极限合成验证');end
    folder=fullfile(out,names{j});writeSU2RocketCase(volume,folder,flight);
    job=embeddedSU2Job('start',folder,struct('pythonExecutable','D:/anacondapython/python.exe','threads',2,'maximumSeconds',60));
    t=tic;while toc(t)<75,[job,receipt]=embeddedSU2Job('poll',job);if job.terminal,break;end;pause(.15);end
    assert(job.terminal && job.supervisorExitCode==0 && receipt.exitCode==0 && ~receipt.timedOut,job.supervisorOutput);receipts{j}=receipt;
    values{j}=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
end
% 预先固定归一化最大误差门槛1e-5；不同编译器舍入不以逐位相同为前提。
errors=zeros(1,3);fields={'pressureRaw','heatFluxRaw','temperatureRaw'};
assert(isequal(values{1}.vertices,values{2}.vertices));
for k=1:3,a=values{1}.(fields{k});b=values{2}.(fields{k});errors(k)=max(abs(a-b))/max(1,max(abs(a)));end
write(fullfile(out,'constant_comparison.json'),struct('maximumNormalizedDifferences',errors,'threshold',1e-5));
assert(all(errors<1e-5),'常量表与官方标量壁不一致。');
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.su2=struct('pythonExecutable','D:/anacondapython/python.exe','iterations',20,'maximumSeconds',60,'outputEvery',20);
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()close(fig));fig.UserData.actions.loadClassic();
fig.UserData.actions.readHybridMesh(meshFolder);state=fig.UserData.workflowState;
assert(isequal(state.mesh.V,surface.V) && isequal(state.mesh.F,surface.F));
assert(contains(fig.UserData.controls.volumeSummary.Text,'124564'));
profile=struct('mode','axialTemperature','x',[0 .158496 .361696 L], ...
    'temperature',[100 150 110 130],'source','v115合成分段线性壁温；非实验数据');
fig.UserData.actions.setWallTemperatureProfile(profile);
assert(strcmp(fig.UserData.controls.thermalMode.Value,'axialTemperature') && strcmp(fig.UserData.controls.wall.Enable,'off'));
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
t=tic;while fig.UserData.workflowState.busy && toc(t)<80,drawnow;pause(.1);end
state=fig.UserData.workflowState;assert(~state.busy && ~isempty(state.researchResult),fig.UserData.status.Text);r=state.researchResult;
folder=fig.UserData.engineJob.folder;applied=r.receipt.thermalProfileApplied;
assert(applied.passed && applied.wallNodes==size(surface.V,1) && ~r.heatFluxIsBoundaryInput);
boundary=resolveSU2ThermalBoundary(surface,profile);assert(isequal(r.wallTemperatureInput,boundary.givenTemperatureByFace));
assert(~r.productionQualified && ~r.thermalBoundary.wallMaterialTemperatureSolved);
ui=fig.UserData.researchUI;ui.select.Value='wallInput';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(.2);waitRender(ui);
assert(ui.geometry.html.UserData.payload.sectionView.valueStatus==2);
assert(isequal(ui.geometry.html.UserData.payload.faceField.values,r.wallTemperatureInput));
ui.inspectFace(100,[]);assert(contains(strjoin(ui.probe.Value,newline),'给定'));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'wall_section.csv'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitRender(ui);
ui.exportSectionTo(fullfile(out,'heat_section.csv'));
assert(max(abs(r.heatFluxRaw-mean(values{2}.heatFluxRaw)))>1); % 实际非零研究场，精确对照由独立审核完成。
report=struct('engineFolder',folder,'constantErrors',errors,'constantThreshold',1e-5, ...
    'constantReceipts',{receipts},'profileApplied',applied,'rendered',ui.geometry.html.UserData.rendered, ...
    'reusedMeshFolder',meshFolder,'noMesherRerun',true,'newQualifiedArea',0,'physicsQualified',false);
save(fullfile(out,'spatial_result.mat'),'r','report','-v7');write(fullfile(out,'report.json'),report);
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('name',fig.Name));fprintf('V115_SPATIAL_UI_READY\n');
    t=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<60,drawnow;pause(.1);end
end
% 更换温度表后旧结果、边界确认必须失效；错误回执表不能附加成新输入。
wrong=boundary;wrong.specification.temperature=wrong.specification.temperature+1;
rejected=false;try,attachSU2ThermalBoundary(r,wrong);catch,rejected=true;end;assert(rejected);
profile.temperature=profile.temperature+1;fig.UserData.actions.setWallTemperatureProfile(profile);
assert(isempty(fig.UserData.workflowState.researchResult) && ~fig.UserData.workflowState.boundaryConfirmed);
write(fullfile(out,'ui_checks.json'),struct('profileEditInvalidatesResult',true,'differentAppliedTableRejected',true));
fprintf('Native constant limit and spatial thermal main workflow passed; no physics qualification.\n');
end
function waitRender(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision&&s.rendered.webglError==0&&~isempty(ui.getSection()));
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
