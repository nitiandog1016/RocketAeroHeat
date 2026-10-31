function lineage=writeSU2CoolPropInitialState(folder,source,python,verifyOnly)
%WRITESU2COOLPROPINITIALSTATE 将同网格定比热源场转换为HEOS初猜，保持原p/T/速度。
% 模块一：新任务逐节点PT/DE校核后写出守恒初态；变物性后局部迭代从零开始。
% 模块二：verifyOnly=true只读复用已完成的转换及SHA，不重复查询、求解或改写历史。
% 模块三：固定使用已测重构保护包；转换不等于守恒迁移、收敛解或物理资格。
if nargin<3,python='python';end
if nargin<4,verifyOnly=false;end
script=fullfile(fileparts(mfilename('fullpath')),'scripts','prepareSU2CoolPropCase.py');
for value={folder,source,python,script}
    assert(~any(ismember(char(value{1}),char([34 10 13]))),'路径包含不支持字符。');
end
p=System.Diagnostics.Process();p.StartInfo.FileName=python;
p.StartInfo.Arguments=sprintf('-X utf8 "%s" --folder "%s" --source "%s" --runtime-flavor rah-coolprop-recon-v3',script,folder,source);
if verifyOnly,p.StartInfo.Arguments=[char(p.StartInfo.Arguments) ' --verify-only'];end
p.StartInfo.UseShellExecute=false;p.StartInfo.CreateNoWindow=true;
p.StartInfo.RedirectStandardOutput=true;p.StartInfo.RedirectStandardError=true;
assert(p.Start());text=char(p.StandardOutput.ReadToEnd());err=char(p.StandardError.ReadToEnd());p.WaitForExit();
assert(p.ExitCode==0,'RocketAeroHeat:SU2:PropertyInitialState','物性初始化被拒绝：%s %s',text,err);
lineage=jsondecode(fileread(fullfile(folder,'coolprop_initialization.json')));
end
