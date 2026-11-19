function check_rejection(base)
%CHECK_REJECTION 核对错误单位、工作区活动边界、跨网格快照和传统工程配置恢复。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
cfg=createCaseConfig();cfg.geometry.stlFile=fullfile(root,'cases','CM011-mm-scale-SI.stl');cfg.geometry.scale=.001;
fig=rocketAeroHeatWorkbench(cfg);cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
fig.Visible='off';press(fig.UserData.controls.load);press(fig.UserData.controls.confirm);
edit(fig.UserData.controls.meshBackend,'cm011Imported');press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;assert(isempty(s.mesh)&&~s.steps(3).available&&~isempty(s.orientedMesh));
assert(contains(fig.UserData.status.Text,'CM-011'),fig.UserData.status.Text);
press(fig.UserData.workflowButtons(1));edit(fig.UserData.controls.scale,1);
press(fig.UserData.controls.load);press(fig.UserData.controls.confirm);
edit(fig.UserData.controls.meshBackend,'cm011Imported');press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;assert(~isempty(s.mesh)&&isfield(s.mesh,'cm011Prepared'),fig.UserData.status.Text);
press(fig.UserData.workflowButtons(3));assert(isempty(evalin('base','flight')));
edit(fig.UserData.controls.cm011Condition,'inclined');press(fig.UserData.controls.confirm);
active=evalin('base','flight');assert(active.alphaDeg==10 && active.speed==.1 && active.wallTemperatureK==300.1);
% 用来自真实前次计算的快照检查完成动作；当前网格的准备回执不同。
s=aeroHeatWorkflowState('beginSolve',fig.UserData.workflowState);
run=strtrim(fileread(fullfile(base,'gui01','run_folder.txt')));prev=load(fullfile(run,'result.mat'),'current');
assert(~strcmp(s.mesh.cm011Prepared.manifestSHA256,prev.current.evidenceIdentity.preparedSHA256));
refused=false;
try
    aeroHeatWorkflowState('completeSolve',s,struct('result',prev.current,'solveId',s.solveId,'sourceRevision',s.revision));
catch ex
    refused=strcmp(ex.identifier,'RAH:CM011:Snapshot');
end
assert(refused);
press(fig.UserData.workflowButtons(2));edit(fig.UserData.controls.meshBackend,'engineering');
assert(isequaln(evalin('base','flight'),cfg.flight));
write=struct('wrongScalePreventedMesh',true,'legacyFlightRestored',true,'activeBoundary',active, ...
    'priorSnapshotRejected',true,'newSolverLaunched',false,'secondMesh',s.mesh.cm011Prepared.directory);
fid=fopen(fullfile(base,'rejections.json'),'w','n','UTF-8');assert(fid>=0);cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(write));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function edit(c,v)
c.Value=v;c.ValueChangedFcn(c,[]);drawnow;
end
