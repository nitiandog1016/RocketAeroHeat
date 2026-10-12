function testBoundaryLayerState()
% 检查零压梯度平板极限，以及强减速时的层流分离提示。
addpath(fileparts(fileparts(mfilename('fullpath'))));
n=101; ds=0.01;
m.F=repmat([1 2 3],n,1);
m.boundaryLayerEdge.method='解析均匀流';
m.boundaryLayerEdge.speed=100*ones(n,1);
m.boundaryLayerEdge.temperature=300*ones(n,1);
m.boundaryLayerEdge.pressure=101325*ones(n,1);
m.boundaryLayerEdge.algebraicallyValid=true(n,1);
p=struct;
p.velocitySource=m.boundaryLayerEdge.method;
p.distance=(0:n-1)'*ds;
p.pathCost=p.distance;
p.predecessorFace=[0;(1:n-1)'];
p.stepLength=[0;ds*ones(n-1,1)];
p.connectedFace=true(n,1);
m.surfaceStreamlines=p;
fs=struct('gamma',1.4,'R',287);
cfg=struct('flight',struct('Twall',300));
uniform=calcBoundaryLayerState(m,fs,cfg);
b=uniform.boundaryLayer;
id=90;
nu=b.referenceKinematicViscosity(id);
Re=100*p.distance(id)/nu;
expected=0.664/sqrt(Re);
assert(abs(b.skinFrictionCoefficient(id)/expected-1)<0.01);
assert(abs(b.pressureGradientParameter(id))<1e-10);
assert(nnz(b.validFace)==n-1);
assert(b.gradientPathSteps==4);
% 一次很小的台阶速度变化不应被单个短面元步长放大为层流分离。
m.boundaryLayerEdge.speed(51:end)=99;
smallStep=calcBoundaryLayerState(m,fs,cfg);
assert(~smallStep.boundaryLayer.possibleLaminarSeparationFace(51));
m.boundaryLayerEdge.speed(51:end)=50;
decelerating=calcBoundaryLayerState(m,fs,cfg);
d=decelerating.boundaryLayer;
assert(d.possibleLaminarSeparationFace(51));
assert(~d.validFace(51) && isnan(d.skinFrictionCoefficient(51)));

% 变比热模式下，Thwaites 参考温度必须与热流模块使用同一焓恢复。
m.boundaryLayerEdge.speed=1500*ones(n,1);
thermal=calcBoundaryLayerState(m,fs,cfg);
thermalB=thermal.boundaryLayer;
expectedTaw=airRecoveryTemperature(300,1500,sqrt(0.72));
assert(strcmp(thermalB.thermalPropertyMode,'thermallyPerfect'));
assert(abs(thermalB.recoveryTemperature(id)-expectedTaw)<1e-8);
assert(abs(thermalB.referenceTemperature(id)- ...
    (300+0.22*(expectedTaw-300)))<1e-8);
assert(contains(thermalB.recoveryTemperatureMethod,'焓恢复'));
caloricCfg=cfg;
caloricCfg.heating.airThermalModel='caloricallyPerfect';
caloric=calcBoundaryLayerState(m,fs,caloricCfg);
caloricB=caloric.boundaryLayer;
caloricTaw=300+sqrt(0.72)*1500^2/ ...
    (2*fs.gamma*fs.R/(fs.gamma-1));
assert(strcmp(caloricB.thermalPropertyMode,'caloricallyPerfect'));
assert(abs(caloricB.recoveryTemperature(id)-caloricTaw)<1e-8);
assert(abs(caloricB.referenceTemperature(id)- ...
    thermalB.referenceTemperature(id))>1);
assert(abs(caloricB.skinFrictionCoefficient(id)- ...
    thermalB.skinFrictionCoefficient(id))>1e-6);
% 模型域外热状态不得用于层流积分；定比热对照模式仍可显式计算。
m.boundaryLayerEdge.temperature=70*ones(n,1);
invalid=calcBoundaryLayerState(m,fs,cfg);
assert(invalid.boundaryLayer.summary.validFaces==0);
assert(invalid.boundaryLayer.summary.thermalPropertyUnavailableFaces==n);
fprintf('Boundary layer passed: flat-plate Cf %.6g (reference %.6g), deceleration lambda %.3f.\n', ...
    b.skinFrictionCoefficient(id),expected,d.pressureGradientParameter(51));
end
