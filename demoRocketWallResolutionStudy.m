function report=demoRocketWallResolutionStudy(folder)
%DEMOROCKETWALLRESOLUTIONSTUDY 重生四档实际近壁分辨率与法向剖面对照工作台。
% 模块用途：从冻结原BE链重放并导出；浏览不隐式启动新求解，不把细网格当资格。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','rocket_wall_resolution_study');end
report=exportRocketWallResolutionStudy(folder);fprintf('完整箭体近壁/剖面对照：%s\n',report.htmlFile);
end
