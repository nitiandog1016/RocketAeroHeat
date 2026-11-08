function launchRocketAeroHeat()
%LAUNCHROCKETAEROHEAT 桌面快捷方式入口，复用正式工作台。
% 模块一：从入口位置解析工程，启动日志按独立MATLAB进程保存。
projectRoot=fileparts(mfilename('fullpath'));
cd(projectRoot);
addpath(projectRoot);
logRoot=fullfile(projectRoot,'outputs','desktop_launcher');
if ~isfolder(logRoot),mkdir(logRoot);end
sessionName=sprintf('session_%s_%d',datestr(now,'yyyymmdd_HHMMSS'),feature('getpid'));
diary(fullfile(logRoot,[sessionName '.log']));
diaryCleanup=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('RocketAeroHeat desktop startup: %s\n',char(datetime('now')));

% 模块二：仅打开工作台，不导入模型或自动执行求解。
main;
appFigure.Icon=fullfile(projectRoot,'assets','icons','rocket-aero-heat.png');
drawnow;
if ~isgraphics(appFigure) || ~strcmp(appFigure.Visible,'on')
    error('RocketAeroHeat:DesktopStartup','工作台窗口未能打开。');
end
state=appFigure.UserData.workflowState;
receipt=struct('status','window-open','timestamp',char(datetime('now')), ...
    'pid',feature('getpid'),'matlabVersion',version,'projectRoot',projectRoot, ...
    'windowName',appFigure.Name,'visible',char(appFigure.Visible), ...
    'workflowButtonCount',numel(appFigure.UserData.workflowButtons), ...
    'currentStep',state.currentStep,'busy',state.busy, ...
    'hasModel',~isempty(state.rawMesh),'hasResult',~isempty(state.result));
receiptPath=fullfile(logRoot,[sessionName '.json']);
fid=fopen(receiptPath,'w','n','UTF-8');
if fid<0
    warning('RocketAeroHeat:DesktopReceipt','无法写入启动记录：%s',receiptPath);
else
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
    clear cleanup;
end
fprintf('RAH_DESKTOP_READY: %s\n',receiptPath);

% 模块三：保留工作台原有关闭/取消回调，关闭窗口后返回快捷方式命令。
% exit由快捷方式调用；在已有MATLAB中手动调用本函数不会退出该会话。
waitfor(appFigure);
fprintf('RAH_DESKTOP_CLOSED\n');
end
