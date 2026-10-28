function report=demoRocketContactStudy(folder)
%DEMOROCKETCONTACTSTUDY 从已确认原BE链重生真实箭体HLLE/HLLC方法对照。
% 同窗切换四档结果/剖面，失败初猜单列；浏览不触发求解，物理资格不自动增加。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','rocket_contact_study');end
report=exportRocketContactStudy(folder);fprintf('完整箭体通量/热流对照：%s\n',report.htmlFile);
end
