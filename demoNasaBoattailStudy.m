function report=demoNasaBoattailStudy(folder)
%DEMONASABOATTAILSTUDY 重算条件压力关联并导出三组名义几何/公开符号科研对照。
% 这是原实验无量纲压力对照，不启动外部CFD、不伪造热流或壁温。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'tests'));
if nargin<1,folder=fullfile(root,'outputs','nasa_boattail_study');end
report=exportNasaBoattailStudy(folder);fprintf('NASA尾锥压力对照：%s\n',report.htmlFile);
end
