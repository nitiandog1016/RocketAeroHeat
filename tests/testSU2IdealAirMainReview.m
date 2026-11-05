function report=testSU2IdealAirMainReview(out)
%TESTSU2IDEALAIRMAINREVIEW 只读兼容及实际窗口观察，复用已完成的三维任务。
% 模块一。旧HEOS/定比热与二维变比热分别按原合同读取。
% 模块二。三维结果从同一主工作台重开，等待实际鼠标检查，不启动求解。
assert(~isfolder(out));mkdir(out);testSU2IdealAirLegacyReads(fullfile(out,'legacy'));
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()close(fig));
for model={'laminar','sst'}
    fig.UserData.actions.readIdealAirReview(fullfile(pwd,'outputs','v154_validation','cases01',model{1}));
    assert(fig.UserData.idealAirReview.summary.dimension==2 && fig.UserData.idealAirReview.convergence.checkpointIdentityConfirmed);
end
clear guard;
source=jsondecode(fileread('outputs/v157_validation/main03/engine_folder.json'));folder=source.folder;
d=load(fullfile(folder,'input_snapshot.mat'),'frozen','budget');
fig=rocketAeroHeatWorkbench(d.frozen.cfg);guard=onCleanup(@()close(fig));
fig.UserData.actions.useSurfaceForSU2(d.frozen.mesh);
fig.UserData.actions.readHybridMesh('outputs/v157_validation/fixture_workbench');
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
ui=fig.UserData.researchUI;ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.setSectionX(.25);ui.inspectFace(2,[]);drawnow;
write(fullfile(out,'ready_for_capture.json'),struct('figureName',fig.Name,'folder',folder,'maximumWaitSeconds',150));
t=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<150,drawnow;pause(.1);end
observed=isfile(fullfile(out,'capture_done.json'));
report=struct('passed',true,'legacyHEOSAndSST',true,'legacyIdealAir2D',true, ...
    'threeDimensionalReopen',true,'nativeMouseObserved',observed,'newNativeRuns',0,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
