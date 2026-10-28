function report=runRocketTemperatureActiveResearch(folder,runFullDomain,initialSourceFolder)
%RUNROCKETTEMPERATUREACTIVERESEARCH 快速校核多温域行，或从显式真实端点重建最多三步。
% 模块一：缺省只运行制造小域双温域选择、冷却行释放和稠密受限子空间校核。
% 模块二：显式true才逐态重建当前NS/C/W/QR及温度行，保留所有原非线性护栏。
% 模块三：输出使用新目录，完整输入随包；未收敛研究不得增加正式求解面积。
% 模块四：可选第三参数选择实际初态目录；旧两参数调用仍缺省Q28，缺失全域输入提前拒绝。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'tests'),fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<1
    base=fullfile(root,'outputs','temperature_active_reproduction');folder=base;index=0;
    while exist(folder,'dir'),index=index+1;folder=sprintf('%s_%02d',base,index);end
end
if nargin<2,runFullDomain=false;end
if nargin<3,initialSourceFolder=fullfile(root,'outputs/v83_validation/krylov28');end
validateattributes(runFullDomain,{'logical'},{'scalar'});
assert((ischar(initialSourceFolder)&&isrow(initialSourceFolder)&&~isempty(initialSourceFolder))|| ...
    (isstring(initialSourceFolder)&&isscalar(initialSourceFolder)&&strlength(initialSourceFolder)>0), ...
    'TemperatureActive:InitialState','初态目录必须为非空路径。');
initialSourceFolder=char(java.io.File(char(initialSourceFolder)).getCanonicalPath());
initialStateFile=fullfile(initialSourceFolder,'wall_velocity_krylov_evidence.mat');
if runFullDomain,assert(isfile(initialStateFile),'TemperatureActive:InitialState','完整研究的真实初态缺失，尚未运行小域或创建输出。');end
assert(~exist(folder,'dir'));mkdir(folder);
report=struct('unit',testRocketTemperatureActiveSet(fullfile(folder,'unit')),'fullDomain',[], ...
    'configuredInitialStateFile',initialStateFile,'configuredInitialStateUsed',runFullDomain, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
if runFullDomain
    report.fullDomain=continueRocketTemperatureActive(fullfile(folder,'continuous'), ...
        initialSourceFolder,3);
end
end
