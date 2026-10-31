function reviewSU2SpatialBoundaryUI(out)
%REVIEWSU2SPATIALBOUNDARYUI 复用本次原结果核对新增边界/网格侧栏，不重建网格或重算。
assert(~isfolder(out));mkdir(out);
report=jsondecode(fileread('outputs/v115_validation/actual01/report.json'));
saved=load(fullfile(report.engineFolder,'input_snapshot.mat'),'frozen');
fig=rocketAeroHeatWorkbench(saved.frozen.cfg);guard=onCleanup(@()close(fig));
fig.UserData.actions.loadClassic();fig.UserData.actions.readHybridMesh(report.reusedMeshFolder);
fig.UserData.actions.setWallTemperatureProfile(saved.frozen.cfg.su2.wallTemperatureProfile);
fig.UserData.actions.confirmBoundary();fig.UserData.actions.readResearch(fullfile(report.engineFolder,'mapped_research.mat'));
assert(~isempty(fig.UserData.workflowState.researchResult));
fig.UserData.actions.selectStep(3);drawnow;
ctl=fig.UserData.controls;assert(strcmp(ctl.thermalMode.Value,'axialTemperature') && strcmp(ctl.wall.Enable,'off'));
fid=fopen(fullfile(out,'ready_for_capture.json'),'w');fprintf(fid,'{"ready":true}');fclose(fid);
fprintf('V115_BOUNDARY_UI_READY\n');t=tic;
while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<60,drawnow;pause(.1);end
fid=fopen(fullfile(out,'checks.json'),'w');fprintf(fid,'{"sameInputReplay":true,"noSolverOrMesher":true,"spatialMode":true,"scalarDisabled":true}');fclose(fid);
end
