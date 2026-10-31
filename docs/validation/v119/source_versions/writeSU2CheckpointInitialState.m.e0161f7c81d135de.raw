function lineage=writeSU2CheckpointInitialState(folder,source,python)
%WRITESU2CHECKPOINTINITIALSTATE 在新任务中复用已确认的同工况原生守恒场。
% 模块一：Python独立核对终态、文件身份、全配置/网格、具名时刻和正状态。
% 模块二：硬链接只读源场，记录本段之前的更新数；伪时间/CFL历史重新开始。
% 模块三：失败在引擎启动前返回，不给不匹配的工况插值或补零。
if nargin<3,python='python';end
script=fullfile(fileparts(mfilename('fullpath')),'scripts','prepareSU2CheckpointContinuation.py');
for item={folder,source,python,script}
    assert(~any(ismember(char(item{1}),char([34 10 13]))),'路径包含不支持字符。');
end
p=System.Diagnostics.Process();p.StartInfo.FileName=python;
p.StartInfo.Arguments=sprintf('-X utf8 "%s" --source "%s" --target "%s"',script,source,folder);
p.StartInfo.UseShellExecute=false;p.StartInfo.CreateNoWindow=true;
p.StartInfo.RedirectStandardOutput=true;p.StartInfo.RedirectStandardError=true;
assert(p.Start());text=char(p.StandardOutput.ReadToEnd());err=char(p.StandardError.ReadToEnd());p.WaitForExit();
assert(p.ExitCode==0,'RocketAeroHeat:SU2:Continuation','检查点准备被拒绝：%s %s',text,err);
lineage=jsondecode(fileread(fullfile(folder,'continuation.json')));
end
