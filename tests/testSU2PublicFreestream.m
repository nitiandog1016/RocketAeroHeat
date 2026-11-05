function report=testSU2PublicFreestream(folder)
%TESTSU2PUBLICFREESTREAM 用MATLAB解析焓熵独立核对公开变比热来流。
% 模块一。保持公开Mach、总温和总压，采用解析h/s核重新求静态。
% 模块二。比较Python准备结果和平均流总态，SST k独立加入守恒能量。
% 模块三。记录局部T/U推导参考，实验总温分母继续等待来源闭合。
target=fullfile(folder,'matlab_freestream.json');assert(~isfile(target));
prepared=jsondecode(fileread(fullfile(folder,'preparation.json')));c=prepared.conditions;
total=idealAirThermodynamicState('PT',c.totalPressure,c.totalTemperature);
temperature=fzero(@(T)energyDifference(T,c.totalPressure,c.mach,total.h), ...
    [60 c.totalTemperature],optimset('TolX',1e-11));
samePressure=idealAirThermodynamicState('PT',c.totalPressure,temperature);
pressure=c.totalPressure*exp((samePressure.s-total.s)/total.gasConstant);
flow=idealAirThermodynamicState('PT',pressure,temperature);
speed=c.mach*flow.sound;
expected=[flow.T flow.p flow.rho flow.h flow.e flow.cp flow.sound speed];
saved=[c.T c.p c.rho c.h c.e c.cp c.sound c.speed];
relativeError=max(abs(expected-saved)./max(1,abs(expected)));
energyError=abs(flow.h+speed^2/2-total.h)/abs(total.h);
entropyError=abs(flow.s-total.s)/max(1,abs(total.s));
assert(relativeError<2e-10 && energyError<2e-11 && entropyError<2e-11);
assert(abs(c.conservativeSpecificEnergy-(flow.e+speed^2/2+c.k))<1e-7);
local=idealAirThermodynamicState('PT',576,81.2);
localTotal=idealAirThermodynamicState('hs',local.h+1274^2/2,local.s);
assert(abs(localTotal.T-prepared.normalization.derivedVariableCpTotalTemperatureK)<1e-7);
report=struct('passed',true,'temperatureK',temperature,'pressurePa',pressure,'speedMS',speed, ...
    'maximumRelativeError',relativeError,'meanTotalEnthalpyRelativeError',energyError, ...
    'entropyRelativeError',entropyError,'derivedLocalTotalTemperatureK',localTotal.T, ...
    'physicsQualified',false);
fid=fopen(target,'w','n','UTF-8');assert(fid>0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
disp(report);
end

function difference=energyDifference(temperature,pressure,mach,enthalpy)
% 有界温区内求平均流总焓方程的根。
value=idealAirThermodynamicState('PT',pressure,temperature);
difference=value.h+.5*(mach*value.sound)^2-enthalpy;
end
