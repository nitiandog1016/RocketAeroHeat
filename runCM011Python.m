function runCM011Python(args,logPath,heartbeat,independent)
%RUNCM011PYTHON 使用独立参数传递路径；计算时可刷新工作台并检查关闭请求。
root=fileparts(mfilename('fullpath'));
if nargin<4,independent=false;end
items=[{'python'},args];argv=javaArray('java.lang.String',numel(items));
for k=1:numel(items),argv(k)=java.lang.String(items{k});end
cmd=java.lang.ProcessBuilder(argv);cmd.directory(java.io.File(root));cmd.redirectErrorStream(true);
env=cmd.environment();env.put('PYTHONIOENCODING','utf-8');
cmd.redirectOutput(java.io.File(logPath));proc=cmd.start();
if nargin>=3 && ~isempty(heartbeat)
    try
        while proc.isAlive(),heartbeat();pause(.2);end
    catch ex
        % 有独立回执的长任务在窗口关闭后仍由自己的监督器完成，可按目录恢复。
        if ~independent,proc.waitFor();end
        rethrow(ex);
    end
end
code=proc.waitFor();
if code~=0
    fid=fopen(logPath,'rt','n','UTF-8');assert(fid>=0,'RAH:CM011:Log','无法读取日志 %s。',logPath);
    closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
    detail=strtrim(fscanf(fid,'%c'));lines=splitlines(string(detail));
    error('RAH:CM011:Process','%s（退出 %d）。日志 %s',char(lines(end)),code,logPath);
end
end
