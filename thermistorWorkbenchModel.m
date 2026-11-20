function [mesh,source]=thermistorWorkbenchModel()
%THERMISTORWORKBENCHMODEL 从校验值已核对的归档取得 B 型珠体的米制顶点和原面。
% 这里只读取几何。工作流在计算步骤另行选择新算或回读的场量。
[folder,clean]=openThermistorArchive(2); %#ok<NASGU>
path=fullfile(folder,'body.json');body=jsondecode(fileread(path));
mesh=struct('V',body.V,'F',body.F);
source=struct('kind','thermistorB','units','m','bodySHA256',cm011FileSHA256(path));
end
