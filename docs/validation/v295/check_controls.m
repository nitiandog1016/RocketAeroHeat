function check_controls(out)
% 模块一：核对 3 m/s 既有热边界与原先的 6 m/s 只读数据包。
[b,fs]=thermistorWorkbenchBoundary();
assert(fs.speed==3 && b.wallTemperatureK==305.15);
[s,r]=runThermistorBeadCase([],'archive',[],[],struct('netHeatingPowerW',0.004));
assert(r.fullSurfaceComputed && abs(r.wholeSurfaceHeatPowerW+.004)<1e-12 && isfield(s,'hybridSolvedWallTemperatureK'));
[frozen,ref]=readThermistorStartup('nearWall6');
assert(ref.fullSurfaceComputed && height(frozen.hybridWall)==5216 && ~isfield(frozen.evidenceIdentity,'runDirectory'));

% 模块二：用短时子进程观察窗口关闭后调用立即返回，子进程仍能独立写出终态。
% 该检查只验收关闭时的进程生命周期；原生求解由五步工作台的新算试验覆盖。
fig=uifigure('Visible','off');guard=onCleanup(@()cleanupFigure(fig)); %#ok<NASGU>
ready=fullfile(out,'child.ready');done=fullfile(out,'child.done');started=tic;closedAt=NaN;
args={'-c',['import pathlib,sys,time,os; p=str(os.getpid()); ' ...
    'pathlib.Path(sys.argv[1]).write_text(p); time.sleep(5); pathlib.Path(sys.argv[2]).write_text(p)'],ready,done};
try
    runCM011Python(args,fullfile(out,'child.log'),@closeAfterStart,true);
    error('RAH:Test:Closed','测试窗口关闭时未返回。');
catch ex
    assert(strcmp(ex.identifier,'RocketAeroHeat:Workbench:Closed'),ex.message);
end
returnedAt=toc(started);assert(isfile(ready) && ~isfile(done) && returnedAt-closedAt<2);
while ~isfile(done) && toc(started)<25,pause(.1);end
assert(isfile(done) && strcmp(fileread(done),fileread(ready)));
p=struct('archiveFaces',height(frozen.hybridWall),'archive56PowerIntoW',ref.wholeSurfaceHeatPowerW, ...
    'steadyPowerIntoW',r.wholeSurfaceHeatPowerW,'childPID',str2double(fileread(done)), ...
    'secondsToReturnAfterClose',returnedAt-closedAt,'childFinishedAfterReturn',true);
fid=fopen(fullfile(out,'controls.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(p));fclose(fid);
    function closeAfterStart()
        if isfile(ready)
            closedAt=toc(started);close(fig);
            error('RocketAeroHeat:Workbench:Closed','工作台已关闭。');
        end
    end
end

function cleanupFigure(fig)
if isgraphics(fig),delete(fig);end
end
