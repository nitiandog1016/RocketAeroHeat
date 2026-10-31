function [job,receipt]=embeddedSU2Job(action,job,options)
%EMBEDDEDSU2JOB 内嵌引擎的异步启动、轮询和取消；不判定物理收敛。
% 模块一：.NET直接启动隐藏监督器，MATLAB/UI不等待长时间数值运算。
% 模块二：取消以唯一令牌提交，由监督器结束它自己拥有的SU2子进程。
% 模块三：只有监督器实际退出且引擎终态存在，才返回terminal=true。
if nargin<3,options=struct;end
receipt=[];
switch lower(action)
    case 'start'
        folder=char(java.io.File(job).getCanonicalPath());
        root=fileparts(mfilename('fullpath'));
        outputRoot=[char(java.io.File(fullfile(root,'outputs')).getCanonicalPath()) filesep];
        assert(startsWith(lower(folder),lower(outputRoot)),'计算目录须在项目outputs内。');
        for name={'case.cfg','mesh.su2'},assert(isfile(fullfile(folder,name{1})));end
        python=option(options,'pythonExecutable','python');
        seconds=option(options,'maximumSeconds',720);threads=option(options,'threads',2);
        validateattributes(seconds,{'numeric'},{'scalar','finite','>=',1,'<=',720});
        validateattributes(threads,{'numeric'},{'scalar','integer','>=',1,'<=',2});
        token=char(java.util.UUID.randomUUID());
        request=struct('workDirectory',folder,'configFile','case.cfg','meshFile','mesh.su2', ...
            'threads',threads,'maximumSeconds',seconds,'dryRun',false,'jobToken',token);
        request.sourceFiles={'readSU2SurfaceResult.m','buildSU2WallProjection.m','mapSU2WallResult.m', ...
            'writeSU2RocketCase.m','embeddedSU2Job.m'};
        if isfile(fullfile(folder,'input_snapshot.mat'))
            request.sourceFiles{end+1}=fullfile(folder,'input_snapshot.mat');
        end
        path=fullfile(folder,'engine_request.json');assert(~isfile(path),'目录已经使用，不能重复启动。');
        write(path,request);
        script=fullfile(root,'scripts','runEmbeddedSU2.py');
        assert(~any(ismember([char(python) script path],char([34 10 13]))));
        process=System.Diagnostics.Process();process.StartInfo.FileName=python;
        process.StartInfo.Arguments=sprintf('-X utf8 "%s" --request "%s"',script,path);
        process.StartInfo.UseShellExecute=false;process.StartInfo.CreateNoWindow=true;
        process.StartInfo.RedirectStandardOutput=true;process.StartInfo.RedirectStandardError=true;
        process.StartInfo.WorkingDirectory=root;
        assert(process.Start(),'无法启动引擎监督器。');
        job=struct('process',process,'folder',folder,'jobToken',token,'terminal',false,'cancelRequested',false);
    case 'poll'
        file=fullfile(job.folder,'engine_process.json');
        if isfile(file),receipt=jsondecode(fileread(file));end
        job.terminal=logical(job.process.HasExited);
        if job.terminal
            job.supervisorExitCode=double(job.process.ExitCode);
            if ~isfield(job,'supervisorOutput')
                job.supervisorOutput=[char(processText(job.process.StandardOutput)) char(processText(job.process.StandardError))];
                fid=fopen(fullfile(job.folder,'supervisor.log'),'w','n','UTF-8');
                if fid>=0,fprintf(fid,'%s',job.supervisorOutput);fclose(fid);end
            end
            if isempty(receipt) || ~isfield(receipt,'status') || ~strcmp(receipt.status,'terminal')
                error('RocketAeroHeat:SU2:Supervisor','监督器已退出但缺少引擎终态，不能接受结果：%s',job.supervisorOutput);
            end
        end
    case 'cancel'
        if ~job.terminal && ~job.process.HasExited
            write(fullfile(job.folder,'engine_cancel.json'),struct('jobToken',job.jobToken));
            job.cancelRequested=true;
        end
    otherwise,error('未知引擎任务动作。');
end
end
function value=processText(stream)
value=stream.ReadToEnd();
end
function value=option(s,key,fallback)
value=fallback;if isfield(s,key),value=s.(key);end
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
