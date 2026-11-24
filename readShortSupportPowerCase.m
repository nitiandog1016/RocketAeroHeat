function [p,archive,seal]=readShortSupportPowerCase(job)
% 模块一：区分版本归档和本次安装生成的量热任务，核对两端功率、目标、候选及实际复算。
root=fullfile(fileparts(mfilename('fullpath')),'cases','short_sphere');
if isempty(job)
    file=fullfile(root,'calorimetry65.json');archive=fullfile(root,'calorimetry65.zip');
else
    file=fullfile(job,'completed.json');archive=fullfile(job,'power65.zip');
end
assert(isfile(file) && isfile(archive),'RAH:Short:Job','所选目录缺少银球量热工况。');
seal=cm011FileSHA256(file);
if isempty(job),assert(strcmp(seal,'94097a483b00142615e272ec6da6baa7f061d106cb9f530f03346c7ff3dcf642'),'RAH:Short:Job','量热归档与发布版本不符。');end
p=jsondecode(fileread(file));b=shortSupportImportBoundary('power',14);
assert(strcmp(p.schema,'short65-sphere-power-v1') && strcmp(p.source.caseSHA256,b.caseSHA256) && ...
    p.faces==32480 && abs(p.areaM2-.000700423853595)<1e-15 && p.experimentalJointAreaM2==0 && ...
    p.targetSphereW==b.targetSphereW && strcmp(cm011FileSHA256(archive),p.archiveSHA256),'RAH:Short:Job','量热任务身份不符。');
lo=p.low.earlyRegionW(1);hi=p.high.earlyRegionW(1);goal=p.targetSphereW;
dt=p.low.excessK+(goal-lo)*(p.high.excessK-p.low.excessK)/(hi-lo);
assert(p.low.excessK==20 && p.high.excessK==40 && lo<goal && goal<hi && abs(hi/2/lo-1)<.002 && ...
    abs(dt-p.balance.excessK)<1e-10 && abs(b.airK+dt-p.temperatureK)<1e-10 && ...
    p.temperatureK==p.balance.temperatureK && abs(p.balance.earlyRegionW(1)/goal-1)<.0002 && ...
    abs(p.balance.lateRegionW(1)/goal-1)<.0002,'RAH:Short:Job','球面功率或独立壁温终点不符。');
assert(strcmp(cm011FileSHA256(file),seal),'RAH:Short:Job','读取期间量热任务回执已改变。');
end
