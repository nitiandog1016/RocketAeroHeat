function report=runRocketTemperatureProjectionResearch(folder,runFullDomain)
%RUNROCKETTEMPERATUREPROJECTIONRESEARCH 快速校核温度投影，或从随包Q26重建真实整箭对照。
% 模块一：缺省执行真实小域温度行、旧内核兼容和稠密受限解校核。
% 模块二：显式runFullDomain=true才重建当前完整系统并对两种128步方向做原NS对照。
% 模块三：缺省输出目录存在时使用新编号；显式目录仍严格拒绝覆盖历史证据。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'tests'),fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<1
    base=fullfile(root,'outputs','temperature_projection_reproduction');folder=base;index=0;
    while exist(folder,'dir'),index=index+1;folder=sprintf('%s_%02d',base,index);end
end
if nargin<2,runFullDomain=false;end
validateattributes(runFullDomain,{'logical'},{'scalar'});
assert(~exist(folder,'dir'),'TemperatureProjection:Overwrite','显式输出目录必须是未存在的新目录。');mkdir(folder);
unit=testRocketTemperatureProjection(fullfile(folder,'unit'));
report=struct('unit',unit,'fullDomainRequested',runFullDomain,'fullDomain',[], ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
if runFullDomain
    state=fullfile(root,'outputs','v83_validation','krylov26','wall_velocity_krylov_evidence.mat');
    geometry=fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat');
    current=fullfile(folder,'current');prepareRocketCurrentPeakSystem(current,state,geometry,'domain');
    report.fullDomain=studyRocketTemperatureProjectedQR(fullfile(folder,'comparison'),fullfile(current,'peak_system.json'),19207);
end
end
