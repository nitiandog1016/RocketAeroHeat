function receipt=runEmbeddedSU2Case(folder,options)
%RUNEMBEDDEDSU2CASE MATLAB同步调用项目内固定SU2，隐藏窗口且限制实际运行时间。
% 模块一：配置和体网格由本项目产生，不要求导入外部CFD结果。
% 模块二：Python只负责进程、预算与文件身份；数值由明确标记的SU2完成。
% 模块三：返回真实终态，不把退出0当作收敛、精度或全面积验证。
if nargin<2,options=struct;end
root=fileparts(mfilename('fullpath'));folder=char(java.io.File(folder).getCanonicalPath());
assert(exist(folder,'dir')==7,'RocketAeroHeat:SU2:Folder','计算目录不存在。');
python='python';if isfield(options,'pythonExecutable'),python=char(options.pythonExecutable);end
threads=2;if isfield(options,'threads'),threads=options.threads;end
seconds=720;if isfield(options,'maximumSeconds'),seconds=options.maximumSeconds;end
validateattributes(threads,{'numeric'},{'scalar','integer','>=',1,'<=',8});
validateattributes(seconds,{'numeric'},{'scalar','finite','positive','<=',720});
request=struct('workDirectory',folder,'configFile','case.cfg','meshFile','mesh.su2', ...
    'threads',threads,'maximumSeconds',seconds,'dryRun',false);
if isfield(options,'dryRun'),request.dryRun=logical(options.dryRun);end
for name={'case.cfg','mesh.su2'}
    assert(exist(fullfile(folder,name{1}),'file')==2,'RocketAeroHeat:SU2:Input','缺少计算输入。');
end
requestFile=fullfile(folder,'engine_request.json');
assert(~exist(requestFile,'file'),'RocketAeroHeat:SU2:Existing','结果目录已经使用，须新建目录。');
fid=fopen(requestFile,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(request));clear guard;
script=fullfile(root,'scripts','runEmbeddedSU2.py');
for text={python,script,requestFile}
    assert(~any(ismember(text{1},char([34 10 13 37]))),'RocketAeroHeat:SU2:Path','执行路径含不支持的字符。');
end
[code,output]=system(sprintf('"%s" "%s" --request "%s"',python,script,requestFile));
receiptFile=fullfile(folder,'engine_process.json');
assert(exist(receiptFile,'file')==2,'RocketAeroHeat:SU2:Launch','未取得引擎回执：%s',output);
receipt=jsondecode(fileread(receiptFile));
assert(strcmp(receipt.status,'terminal'),'RocketAeroHeat:SU2:Running','引擎尚未取得终态。');
assert(code==0 && receipt.exitCode==0 && ~receipt.timedOut, ...
    'RocketAeroHeat:SU2:Failed','SU2未正常完成，详见%s；%s',receiptFile,output);
fprintf('Embedded SU2 terminal exit0, %.3f s; physics qualification is separate.\n',receipt.elapsedSeconds);
end
