function report=runRocketAugmentedQRResearch(folder,runFullDomain)
%RUNROCKETAUGMENTEDQRRESEARCH 自包含运行分块QR校核，或从随包真实Q25重建整箭研究。
% 模块一：缺省只运行小域全部列/转置/零与非零初猜稠密校核，适合快速检查安装。
% 模块二：显式runFullDomain=true才重建真实217272CV当前系统并进行原NS门槛下的一步试算。
% 模块三：工程包保留完整真实输入；不会覆盖历史结果或继承旧状态的方向资格。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'tests'),fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<1,folder=fullfile(root,'outputs','augmented_qr_reproduction');end
if nargin<2,runFullDomain=false;end
validateattributes(runFullDomain,{'logical'},{'scalar'});
assert(~exist(folder,'dir'),'AugmentedQR:Overwrite','请使用未存在的新输出目录。');mkdir(folder);
unit=testRocketAugmentedPatchQR(fullfile(folder,'unit'));
report=struct('unit',unit,'fullDomainRequested',runFullDomain,'fullDomain',[], ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
if runFullDomain
    state=fullfile(root,'outputs','v83_validation','krylov25','wall_velocity_krylov_evidence.mat');
    geometry=fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat');
    assert(isfile(state)&&isfile(geometry),'AugmentedQR:Input','工程包缺少真实Q25或原网格。');
    current=fullfile(folder,'current');prepareRocketCurrentPeakSystem(current,state,geometry,'domain');
    report.fullDomain=studyRocketAugmentedPatchQR(fullfile(folder,'full'),fullfile(current,'peak_system.json'));
end
end
