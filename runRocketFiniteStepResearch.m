function report=runRocketFiniteStepResearch(folder)
%RUNROCKETFINITESTEPRESEARCH 随工程包交付的自包含小复现入口，不需要旧机缓存路径。
% 保留真实Q22/Q24几何和状态；运行六步同预算比较及独立代数守门。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'tests'),fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<1,folder=fullfile(root,'outputs',['finite_step_replay_' datestr(now,'yyyymmdd_HHMMSS')]);end
assert(~exist(folder,'dir'),'FiniteStep:Overwrite','请使用新目录。');mkdir(folder);
fixture=fullfile(root,'outputs/v84_validation/diagnosis01/finite_step_patch.mat');assert(isfile(fixture),'FiniteStep:Fixture','交付夹具缺失。');
unit=testRocketMinimumStateDirection(fullfile(folder,'unit'));
comparison=studyRocketFiniteStepPatch(fullfile(folder,'comparison'),fixture,'minimum-state');
matched=studyRocketFiniteStepPatch(fullfile(folder,'matched'),fixture,'minimum-state-matched');
report=struct('status','completed-packaged-real-geometry-local-reproduction','unit',unit,'comparison',comparison,'matchedConstraintComparison',matched, ...
    'wholeRocketRecomputed',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'research_replay.mat'),'report','-v7');
end
