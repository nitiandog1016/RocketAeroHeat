function report=demoRocketAxisymmetricSensitivity(folder)
%DEMOROCKETAXISYMMETRICSENSITIVITY 从已交付证据重生四档定常箭体对照工作台。
% 仅逐链检查与导出，不隐式启动长时求解；返回htmlFile供离线打开。
% 当前热流网格敏感性明显，数值定常不等于网格独立或生产资格。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','rocket_axisymmetric_sensitivity');end
report=exportRocketAxisymmetricSensitivity(folder);
fprintf('完整箭体网格/外域对照：%s\n',report.htmlFile);
end
