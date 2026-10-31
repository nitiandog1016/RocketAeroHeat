function report=auditSU2PropertyAssumption(out)
%AUDITSU2PROPERTYASSUMPTION 量化当前定比热假设与已有NASA7空气模块的差异。
% 这是物性模型对照，不是流场复算或实验验证；不能用后处理改温度替代能量方程。
assert(~isfile(out));s=load('outputs/v118_validation/actual01/surface.mat','fs');fs=s.fs;
R=287.05287;gamma=1.4;constantCp=gamma*R/(gamma-1);
T=[fs.temperature;300;600;1000;1359.0545016679803;1619.8676846448952];
[cp,h]=airThermalProperties(T);[~,hinf]=airThermalProperties(fs.temperature);
target=hinf+fs.speed^2/2;
T0=fzero(@(t)enthalpy(t)-target,[fs.temperature 2500]);
T0perfect=fs.temperature+fs.speed^2/(2*constantCp);
report=struct('comparisonOnly',true,'physicsQualified',false,'temperatureK',T, ...
    'constantCpJkgK',constantCp,'nasa7CpJkgK',cp,'nasa7CpOverConstantMinus1',cp/constantCp-1, ...
    'nasa7EnthalpyRiseFromFreestreamJkg',h-hinf, ...
    'constantCpEnthalpyRiseFromFreestreamJkg',constantCp*(T-fs.temperature), ...
    'sameVelocityAndTinf',true,'caloricallyPerfectTotalTemperatureK',T0perfect, ...
    'nasa7EnthalpyTotalTemperatureK',T0,'totalTemperatureRelativeDifference',T0perfect/T0-1, ...
    'note','NASA7既有固定组分对照；不是本次SU2已使用的物性，未重算p/q或授予资格。');
fid=fopen(out,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Property-model comparison: constant-cp T0 %.6f K, existing NASA7 T0 %.6f K.\n',T0perfect,T0);
end
function h=enthalpy(T)
[~,h]=airThermalProperties(T);
end
