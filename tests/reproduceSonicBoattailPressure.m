function report=reproduceSonicBoattailPressure(folder)
%REPRODUCESONICBOATTAILPRESSURE 保存声速流线关联复算，不以此冒充实验验证。
% 模块一：由同一给定柱压计算原报告角度的条件候选，并保留全部未知热量。
% 模块二：拒绝非法角度/柱压，另由Python积分形式独立核对数学结果。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(nargin==1&&~exist(folder,'dir'),'SonicBoattailStudy:Overwrite','须指定新的复算目录。');
angles=0:3:21;candidate=estimateSonicBoattailPressure(1000,angles);
assert(candidate.candidatePressureRatio(1)==1&&candidate.sonicStreamlineMachProxy(1)==1);
assert(all(diff(candidate.candidatePressureRatio)<0)&&all(diff(candidate.sonicStreamlineMachProxy)>0));
assert(all(isnan(candidate.heatFluxWm2))&&all(isnan(candidate.edgeTemperatureK))&& ...
    all(isnan(candidate.recoveryTemperatureK))&&all(isnan(candidate.wallTemperatureSolvedK))&& ...
    ~candidate.productionQualified&&~candidate.externalEdgeStateSolved);
bad={@()estimateSonicBoattailPressure(0,6),@()estimateSonicBoattailPressure(NaN,6), ...
    @()estimateSonicBoattailPressure(1000,-1),@()estimateSonicBoattailPressure(1000,24), ...
    @()estimateSonicBoattailPressure([1 2],6)};
for k=1:numel(bad)
    rejected=false;try,bad{k}();catch,rejected=true;end
    assert(rejected,'SonicBoattailStudy:Input','非法工况未被拒绝。');
end
report=struct('status','mathematical-correlation-reproduced-awaiting-independent-quadrature', ...
    'candidate',candidate,'inputRejections',numel(bad),'experimentalComparisonPerformed',false, ...
    'physicalValidationComplete',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'sonic_pressure_evidence.mat'),'report');
fid=fopen(fullfile(folder,'sonic_pressure_reproduction.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('声速流线压力关联已复算，p6/p_c=%.8f，p12/p_c=%.8f；不是实验验证。\n', ...
    candidate.candidatePressureRatio(3),candidate.candidatePressureRatio(5));
end
