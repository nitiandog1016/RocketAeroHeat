function report=demoRocketExactFluxStudy(folder,evidenceRoot,resultFolder)
%DEMOROCKETEXACTFLUXSTUDY 重建真实原面通量/候选路径对照，不自动启动长时计算。
% 模块一：三个状态保持同一真实几何、工况和源初猜，分别保留离散方程身份。
% 模块二：可显式选择保守或primitive研究端点；未知温度和未收敛状态保持。
root=fileparts(mfilename('fullpath'));
if nargin<1,folder=fullfile(root,'outputs','rocket_exact_flux_study');end
if nargin<2
    evidenceRoot=fullfile(root,'outputs','v82_validation');
    if ~exist(fullfile(evidenceRoot,'integrate01','exact_integration_evidence.mat'),'file'),evidenceRoot=fullfile(root,'validation','v82');end
end
if nargin<3,resultFolder='exact_krylov01';end
report=exportRocketExactFluxStudy(folder,evidenceRoot,resultFolder);
fprintf('真实20度精确通量/三维NS研究：%s\n',report.htmlFile);
end
