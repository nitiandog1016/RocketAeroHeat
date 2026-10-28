function report=demoRocketShellStudy(folder,evidenceRoot)
%DEMOROCKETSHELLSTUDY 导出已保存的真实20度三维启动研究，不自动长跑或冒充收敛。
% 模块用途：同窗比较旧失败初态、Dirichlet修正初态和真实三步；无外部CFD输入。
root=fileparts(mfilename('fullpath'));
if nargin<1,folder=fullfile(root,'outputs','rocket_shell_study');end
if nargin<2
    evidenceRoot=fullfile(root,'outputs','v80_validation');
    if ~exist(fullfile(evidenceRoot,'advance01','advance_evidence.mat'),'file')
        evidenceRoot=fullfile(root,'validation','v80');
    end
end
report=exportRocketShellStudy(folder,evidenceRoot);fprintf('真实20度三维数值研究：%s\n',report.htmlFile);
end
