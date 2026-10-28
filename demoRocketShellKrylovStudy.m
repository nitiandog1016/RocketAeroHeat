function report=demoRocketShellKrylovStudy(folder,evidenceRoot)
%DEMOROCKETSHELLKRYLOVSTUDY 导出真实20度五组功策略/迭代结果，不自动长跑或假定稳态。
root=fileparts(mfilename('fullpath'));
if nargin<1,folder=fullfile(root,'outputs','rocket_shell_krylov_study');end
if nargin<2
    evidenceRoot=fullfile(root,'outputs','v81_validation');
    if ~exist(fullfile(evidenceRoot,'arithmetic01','steady_krylov_evidence.mat'),'file'),evidenceRoot=fullfile(root,'validation','v81');end
end
report=exportRocketShellKrylovStudy(folder,evidenceRoot);fprintf('真实20度功速度/Krylov研究：%s\n',report.htmlFile);
end
