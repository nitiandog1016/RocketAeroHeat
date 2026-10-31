function [cp,h,s0]=shockLayerThermalBatch(T,gas)
%SHOCKLAYERTHERMALBATCH 同一激波层物性约定的批量cp、焓和参考压力下熵。
% 模块一：严格检查声明温域，既不外推也不裁温。
% 模块二：定比热保留cp*T；NASA7使用项目1000K连续积分约定。
% 模块三：s0对应gas.referencePressurePa，实际压力熵另减R*log(p/pref)。
assert(isreal(T)&&all(isfinite(T),'all')&& ...
    all(T>=gas.temperatureRangeK(1)&T<=gas.temperatureRangeK(2),'all'), ...
    'ShockPatch:Temperature','激波层批量物性超出声明温域，禁止裁温或外推。');
if strcmp(gas.thermalPropertyMode,'caloricallyPerfect')
    cp=repmat(gas.gamma*gas.R/(gas.gamma-1),size(T));h=cp.*T;
    s0=cp.*log(T/gas.referenceTemperatureK);
else
    assert(strcmp(gas.thermalPropertyMode,'thermallyPerfect')&& ...
        strcmp(gas.thermalPropertyConvention,'continuous-project-cp-v1')&& ...
        abs(gas.R-287.05)<=5&&gas.temperatureRangeK(1)>=80&&gas.temperatureRangeK(2)<=2500, ...
        'ShockPatch:Gas','NASA7须使用项目固定组分空气的连续积分约定。');
    [cp,h,s0]=surfaceFlowAirThermalProperties(T);
    [~,~,reference]=surfaceFlowAirThermalProperties(gas.referenceTemperatureK);
    s0=s0-reference;
end
assert(all(isfinite(cp)&cp>gas.R&isfinite(h)&isfinite(s0),'all'), ...
    'ShockPatch:Gas','批量物性失效，不能切换至另一种气体模型。');
end
