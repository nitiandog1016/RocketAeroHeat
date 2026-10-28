%% RocketAeroHeat 科研工作台入口
% 模块一：准备工程路径和配置；启动时不导入模型、不猜测朝向或计算场量。
% 模块二：五步骤侧边栏在同一窗口管理模型/网格/边界/显式计算/结果。
% 原逐段脚本完整保留为 mainLegacy.m，供既有流程对照，不是默认入口。
projectRoot=fileparts(mfilename('fullpath'));
addpath(projectRoot);
cfg=createCaseConfig();
appFigure=rocketAeroHeatWorkbench(cfg);
