function report=runRocketTemperatureLookaheadResearch(folder,runFullDomain,initialSourceFolder)
%RUNROCKETTEMPERATURELOOKAHEADRESEARCH 复现前瞻选面小域，或从显式真实端点重建一次研究方向。
% 模块一：缺省运行随包上下温域制造小域，复算原整带与目标驱动集合，0新NS。
% 模块二：显式true才重建真实当前NS/C/W/QR与自由方向，原温域/NS/下降门槛不变。
% 模块三：这是未收敛研究入口；全面积物理资格、网格/外域/公开验证分别验收。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'tests'),fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<1
    base=fullfile(root,'outputs','temperature_lookahead_reproduction');folder=base;index=0;
    while exist(folder,'dir'),index=index+1;folder=sprintf('%s_%02d',base,index);end
end
if nargin<2,runFullDomain=false;end
if nargin<3,initialSourceFolder=fullfile(root,'outputs/v83_validation/krylov38');end
validateattributes(runFullDomain,{'logical'},{'scalar'});
initialSourceFolder=char(java.io.File(char(initialSourceFolder)).getCanonicalPath());
state=fullfile(initialSourceFolder,'wall_velocity_krylov_evidence.mat');
if runFullDomain,assert(isfile(state),'TemperatureLookahead:InitialState','真实初态缺失，尚未创建输出或执行小域。');end
assert(~exist(folder,'dir'));mkdir(folder);
report=struct('unit',testRocketTemperatureTargetedSmall(fullfile(folder,'unit')), ...
    'fullDomain',[],'configuredInitialStateFile',state,'configuredInitialStateUsed',runFullDomain, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
if runFullDomain
    report.fullDomain=studyRocketTemperatureLookahead(fullfile(folder,'continuous'),initialSourceFolder,1/8, ...
        'original-and-predicted-required-faces');
end
end
