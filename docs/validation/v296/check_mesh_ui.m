function check_mesh_ui(out,warm)
% 模块一：在一个没有既有网格或初场的安装中，从五步侧栏开始全周网格准备。
fig=rocketAeroHeatWorkbench();fig.Visible='off';cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.beadModel.ButtonPushedFcn(c.beadModel,[]);
c=fig.UserData.controls;c.beadCells.Value=48;c.beadCells.ValueChangedFcn(c.beadCells,[]);
root=fileparts(which('rocketAeroHeatWorkbench'));cache=fullfile(root,'outputs','thermistor_inputs','B6-56');
if nargin<2,warm=false;end
assert(isfolder(cache)==warm);pulses=0;
if warm,before=dir(fullfile(cache,'mesh48-*','debugger_process.json'));end
% 只有独立的几何进程尚在工作时才计数。从 UI 事件循环取样，核实首次建网期间可响应。
watch=timer('ExecutionMode','fixedSpacing','Period',1,'TimerFcn',@sample);
watchCleanup=onCleanup(@()disposeTimer(watch)); %#ok<NASGU>
start(watch);c=fig.UserData.controls;c.generate.ButtonPushedFcn(c.generate,[]);stop(watch);
clear watchCleanup;
s=fig.UserData.workflowState;assert(~isempty(s.mesh),fig.UserData.status.Text);p=s.mesh.thermistorPrepared;
assert(p.nodes==172482 && p.cells==172032 && size(p.preview.F,1)==6048);
if warm
    assert(pulses==0 && isequal(before,dir(fullfile(cache,'mesh48-*','debugger_process.json'))));
else
    assert(pulses>0);
end
assert(~isfile(fullfile(p.directory,'seed_00000.csv')) && isempty(s.result) && ~s.steps(4).available);
proof=struct('warmMeshCache',warm,'eventCallbacksDuringGridWork',pulses,'nodes',p.nodes,'cells',p.cells, ...
    'previewTriangles',size(p.preview.F,1),'nativeSeedCreated',false,'meshSHA256',p.meshSHA256);
fid=fopen(fullfile(out,'mesh_ui.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
    function sample(~,~)
        logs=dir(fullfile(cache,'mesh48-*','debugger_process.json'));
        if isempty(logs),return;end
        r=jsondecode(fileread(fullfile(logs(1).folder,logs(1).name)));
        if strcmp(r.status,'running'),pulses=pulses+1;end
    end
end

function disposeTimer(watch)
if isvalid(watch),stop(watch);delete(watch);end
end
