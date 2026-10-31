function gas=createSurfaceFlowGas(mode,R,gamma,temperatureRangeK)
%CREATESURFACEFLOWGAS 明确创建曲面流线使用的热物性，不进行静默回退。
% 模块一：定比热沿用原gamma/R；变比热调用项目固定组分干空气NASA7。
% 模块二：冻结温域及参考温压，保留焓原零点和热物性来源。
if nargin<4
    if strcmp(mode,'thermallyPerfect'),temperatureRangeK=[80 2500];else,temperatureRangeK=[80 2000];end
end
assert(any(strcmp(mode,{'caloricallyPerfect','thermallyPerfect'})),'SurfaceFlow:Gas','未知热物性模式。');
validateattributes(R,{'double'},{'scalar','real','finite','positive'});
validateattributes(gamma,{'double'},{'scalar','real','finite','>',1});
validateattributes(temperatureRangeK,{'double'},{'real','finite','size',[1 2],'positive'});
assert(temperatureRangeK(2)>temperatureRangeK(1),'SurfaceFlow:Gas','声明温域须严格递增。');
gas=struct('thermalPropertyMode',mode,'gamma',gamma,'R',R, ...
    'temperatureRangeK',temperatureRangeK,'referenceTemperatureK',300,'referencePressurePa',101325);
if strcmp(mode,'thermallyPerfect')
    gas.thermalPropertyConvention='continuous-project-cp-v1';
    gas.source='项目固定组分NASA7/低温冻结cp；1000K处h/s连续匹配积分约定';
else
    gas.thermalPropertyConvention='calorically-perfect-v94';
    gas.source='明确给定gamma/R的定比热完全气体；不是NASA7失效回退';
end
surfaceFlowThermodynamicState(gas.referenceTemperatureK,gas.referencePressurePa,gas);
end
