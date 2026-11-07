function result = runEilmerCase(action, options)
% runEilmerCase  调用固定 Eilmer 后端；准备、运行、读取必须分别明确调用。
% 本文件是完整 MATLAB 入口。不会自动安装环境、接续旧场或取得物理资格。
% action: input / prepare / development / run / read / inspect；options 字段见各模块。
% 公共字段：caseDirectory、recordDirectory；seconds 默认60，outputMiB默认128。
% input另需 meshFile、condition；run另需 preparedReceipt；read另需 loadsFile、runReceipt。
% development另需sourceCase、sourcePreparedReceipt、sourceGeometryAudit、condition。
% inspect另需runReceipt；只读返回失败状态，不调用求解器，不把exit0当物理/过程完成。

%% 模块一：定位项目和固定运行时，拒绝覆盖已有回执。
arguments
    action (1,1) string
    options (1,1) struct
end
root = fileparts(mfilename('fullpath'));
assert(ispc, 'RocketAeroHeat:EilmerHost', '此入口面向已安装的 Windows/WSL 运行时。');
assert(any(action == ["input","prepare","development","run","read","inspect"]), ...
    'RocketAeroHeat:EilmerAction', '未知后端操作。');
caseDir = absolutePath(required(options,'caseDirectory'));
recordDir = absolutePath(required(options,'recordDirectory'));
assert(~isfolder(recordDir) && ~isfile(recordDir), ...
    'RocketAeroHeat:EilmerOverwrite', '回执目录已存在，请使用新目录。');
seconds = setting(options,'seconds',60);
outputMiB = setting(options,'outputMiB',128);
validateattributes(seconds,{'numeric'},{'scalar','real','finite','positive','<=',1800});
validateattributes(outputMiB,{'numeric'},{'scalar','integer','positive','<=',2048});
mkdir(recordDir);

%% 模块二：构造参数数组；不拼接 shell，不在浏览或读取时暗中求解。
switch action
    case "input"
        condition = required(options,'condition');
        conditionFile = fullfile(recordDir,'condition.json');
        writeJSON(conditionFile,condition);
        script = 'prepareEilmerRocketCase.py';
        args = {'--mesh',wslPath(absolutePath(required(options,'meshFile'))), ...
            '--output',wslPath(caseDir),'--condition-json',wslPath(conditionFile)};
        reportFile = fullfile(caseDir,'case_manifest.json');
    case {"prepare","run"}
        script = 'runPreparedEilmerCase.py';
        args = {'--case',wslPath(caseDir),'--stage',char(action), ...
            '--record',wslPath(fullfile(recordDir,'native')), ...
            '--seconds',num2str(seconds,17),'--output-mib',num2str(outputMiB)};
        if action == "run"
            args = [args, {'--prepared-receipt', ...
                wslPath(absolutePath(required(options,'preparedReceipt')))}];
        end
        reportFile = fullfile(recordDir,'native','report.json');
    case "development"
        % 复用已审核原生网格；新初场由原生写入并与原均匀初态逐字核对。
        conditionFile=fullfile(recordDir,'condition.json');
        writeJSON(conditionFile,required(options,'condition'));
        script='prepareEilmerDevelopmentCase.py';
        args={'--source',wslPath(absolutePath(required(options,'sourceCase'))), ...
            '--prepared',wslPath(absolutePath(required(options,'sourcePreparedReceipt'))), ...
            '--audit',wslPath(absolutePath(required(options,'sourceGeometryAudit'))), ...
            '--case',wslPath(caseDir),'--record',wslPath(fullfile(recordDir,'native')), ...
            '--condition',wslPath(conditionFile),'--seconds',num2str(seconds,17), ...
            '--output-mib',num2str(outputMiB)};
        reportFile=fullfile(recordDir,'native','report.json');
    case "inspect"
        % 工作台可显示明确失败原因；没有原生终态时禁止伪装成空的成功云图。
        script='inspectEilmerRun.py';reportFile=fullfile(recordDir,'inspection.json');
        args={'--case',wslPath(caseDir),'--receipt',wslPath(absolutePath(required(options,'runReceipt'))), ...
            '--out',wslPath(reportFile)};
    case "read"
        script = 'readEilmerNative.py';
        reportFile = fullfile(recordDir,'read-report.json');
        args = {'--case',wslPath(caseDir),'--report',wslPath(reportFile), ...
            '--loads',wslPath(absolutePath(required(options,'loadsFile'))), ...
            '--run-receipt',wslPath(absolutePath(required(options,'runReceipt'))), ...
            '--csv',wslPath(fullfile(recordDir,'wall.csv'))};
end
wsl = fullfile(root,'tools','runtime','windows-wsl','wsl.exe');
assert(isfile(wsl), 'RocketAeroHeat:EilmerRuntime', '固定 WSL 运行时不存在。');
argv = [{wsl,'-d','RocketAeroHeat-Ubuntu-24.04','-u','rah','--exec', ...
    '/usr/bin/timeout','--signal=TERM','--kill-after=5',num2str(ceil(seconds)+10), ...
    'python3',wslPath(fullfile(root,'scripts',script))},args];

%% 模块三：无窗口等待；Linux监督负责进程组/时间/容量和输入身份。
% MATLAB取消时关闭本次Windows入口；Linux独立timeout仍限制已启动过程寿命。
list = java.util.ArrayList();
for i = 1:numel(argv), list.add(java.lang.String(argv{i})); end
builder = java.lang.ProcessBuilder(list);
builder.directory(java.io.File(root));
builder.redirectErrorStream(true);
builder.redirectOutput(java.io.File(fullfile(recordDir,'launcher.log')));
timer = tic;
process = builder.start();
cleanup = onCleanup(@() closeProcess(process));
while process.isAlive()
    pause(0.1); drawnow limitrate;
    if toc(timer) > seconds + 35
        process.destroyForcibly();
        error('RocketAeroHeat:EilmerTimeout','WSL入口超时；详见本次回执目录。');
    end
end
code = process.exitValue();
launcher = struct('action',action,'exitCode',double(code),'elapsedSeconds',toc(timer), ...
    'reportFile',reportFile,'physicsQualified',false,'newQualifiedAreaM2',0);
writeJSON(fullfile(recordDir,'launcher.json'),launcher);
assert(code == 0 && isfile(reportFile), 'RocketAeroHeat:EilmerFailed', ...
    '后端未通过；保留日志：%s', fullfile(recordDir,'launcher.log'));
result = jsondecode(fileread(reportFile));
if isfield(result,'passed') && action ~= "inspect"
    assert(result.passed,'RocketAeroHeat:EilmerRejected','后端拒绝结果。');
end
result.matlabLauncher = launcher;
% 此处不写工程结果状态、不将研究CSV接为已验证全面积解。
end

function value = required(options,key)
assert(isfield(options,key),'RocketAeroHeat:EilmerOption','缺少字段：%s',key);
value = options.(key);
end

function value = setting(options,key,fallback)
if isfield(options,key), value=options.(key); else, value=fallback; end
end

function path = absolutePath(value)
path = char(string(value));
assert(isrow(path) && ~isempty(regexp(path,'^[A-Za-z]:[\\/]','once')) && ...
    ~any(double(path)<32),'RocketAeroHeat:EilmerPath','须提供本地盘绝对路径。');
path = char(java.io.File(path).getCanonicalPath());
end

function path = wslPath(value)
path = ['/mnt/',lower(value(1)),'/',strrep(value(4:end),'\','/')];
end

function writeJSON(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'无法写入回执。');
closer=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value));
end

function closeProcess(process)
if process.isAlive(), process.destroy(); end
end
