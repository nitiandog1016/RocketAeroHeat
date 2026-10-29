function [rhs,info]=pressureSurfaceFlowRHS(state,problem,options)
%PRESSURESURFACEFLOWRHS 在一般光滑曲面坐标上积分切向Euler流向方程。
% 模块一：读取真实曲面的一阶/二阶导数，构造尺度无关的切向投影与联络。
% 模块二：核对同一理想气体的p/rho/T/U及导数、沿程总焓与熵；输入不补值。
% 模块三：压力横向梯度决定转向，曲面联络负责几何输运；不求压力或热流。
% 理论：D_s t=-(I-t*t')*grad_surface(p)/(rho*U^2)。这是物面Euler
% 横向动量关系；参见Riley/DeJarnette NTRS 20040121205印刷第4页式13。
% 本函数不是该文完整激波层/VII算法，也不以给定压力代替压力求解。
if nargin<3, options=pressureSurfaceFlowOptions; end
validateattributes(state,{'double'},{'real','finite','size',[4 1]});
sample=problem.sample(state(1:2));
assert(isstruct(sample)&&isscalar(sample),'SurfaceFlow:Sample','采样器须返回标量结构。');
assert(isfield(sample,'isAttached')&&isequal(sample.isAttached,true)&& ...
    isfield(sample,'isSmooth')&&isequal(sample.isSmooth,true), ...
    'SurfaceFlow:Unsupported','采样点未声明处于同一光滑附着区；不跨越分离或压力接口。');
required={'position','firstDerivatives','secondDerivatives','pressure','density','temperature','speed', ...
    'pressureGradient','densityGradient','temperatureGradient','speedGradient'};
assert(all(isfield(sample,required)),'SurfaceFlow:Sample','缺少曲面、物性或两向导数字段。');
A=sample.firstDerivatives;B=sample.secondDerivatives;
validateattributes(sample.position,{'double'},{'real','finite','size',[3 1]});
validateattributes(A,{'double'},{'real','finite','size',[3 2]});
validateattributes(B,{'double'},{'real','finite','size',[3 3]});
for key={'pressure','density','temperature'}
    validateattributes(sample.(key{1}),{'double'},{'real','finite','scalar','positive'});
end
validateattributes(sample.speed,{'double'},{'real','finite','scalar','nonnegative'});
for key={'pressureGradient','densityGradient','temperatureGradient','speedGradient'}
    validateattributes(sample.(key{1}),{'double'},{'real','finite','size',[2 1]});
end
g=problem.gas;R=g.R;
thermal=surfaceFlowThermodynamicState(sample.temperature,sample.pressure,g);cp=thermal.specificHeat;
assert(sample.temperature>=g.temperatureRangeK(1)&&sample.temperature<=g.temperatureRangeK(2), ...
    'SurfaceFlow:Temperature','输入外缘温度超出声明物性范围；不裁温。');
assert(sample.speed>options.minimumSpeedMS,'SurfaceFlow:Stagnation','驻点奇点须另行启动；不猜速度。');

%% 一、缩放曲面基，避免把坐标单位差异误判为几何退化
lengths=vecnorm(A,2,1)';assert(all(lengths>0),'SurfaceFlow:Chart','曲面基退化。');
unitA=A./lengths';C=unitA'*unitA;normal=cross(unitA(:,1),unitA(:,2));
sine=norm(normal);assert(sine>=options.minimumChartAngleSin,'SurfaceFlow:Chart','曲面坐标接近退化。');
normal=normal/sine;v=state(3:4);physical=A*v;normT=norm(physical);
assert(isfinite(normT)&&normT>0,'SurfaceFlow:Direction','流线切向量为零或非有限。');
tangent=physical/normT;coordinateUnit=v/normT;
physicalGradient=@(covector)unitA*(C\(covector./lengths));
gradP=physicalGradient(sample.pressureGradient);gradRho=physicalGradient(sample.densityGradient);
gradT=physicalGradient(sample.temperatureGradient);gradU=physicalGradient(sample.speedGradient);

%% 二、输入状态及沿程Euler/能量/熵合同；不是外流或分离验证
p=sample.pressure;rho=sample.density;T=sample.temperature;U=sample.speed;
assert(abs(p-rho*R*T)<=options.relativeStateTolerance*max(p,rho*R*T), ...
    'SurfaceFlow:EOS','输入p/rho/T不属于同一状态方程。');
eosDerivative=gradP-R*(T*gradRho+rho*gradT);
eosScale=max([norm(gradP),R*T*norm(gradRho),R*rho*norm(gradT)]);
assert(norm(eosDerivative)<=options.absolutePressureGradientTolerancePaPerM+ ...
    options.relativeGradientTolerance*eosScale,'SurfaceFlow:EOSGradient','输入状态与两向导数不一致。');
momentum=tangent'*(rho*U*gradU+gradP);
momentumScale=max(norm(gradP),rho*U*norm(gradU));
assert(abs(momentum)<=options.absolutePressureGradientTolerancePaPerM+ ...
    options.relativeGradientTolerance*momentumScale,'SurfaceFlow:Momentum','沿流线速度与压力梯度不相容。');
energyGradient=cp*gradT+U*gradU;
entropyGradient=cp*gradT/T-R*gradP/p;
assert(abs(tangent'*energyGradient)<=options.absoluteEnergyGradientToleranceJPerKgM+ ...
    options.relativeGradientTolerance*max(cp*norm(gradT),U*norm(gradU)), ...
    'SurfaceFlow:EnergyGradient','输入沿程总焓导数不闭合。');
assert(abs(tangent'*entropyGradient)<=options.absoluteEntropyGradientToleranceJPerKgKM+ ...
    options.relativeGradientTolerance*max(cp*norm(gradT)/T,R*norm(gradP)/p), ...
    'SurfaceFlow:EntropyGradient','输入沿程熵不守恒；不能逐站重置总压。');
totalH=thermal.enthalpy+0.5*U^2;entropy=thermal.entropy;
if isfield(problem,'referenceTotalEnthalpy')
    assert(abs(totalH-problem.referenceTotalEnthalpy)<=options.relativeStateTolerance* ...
        max(1,abs(problem.referenceTotalEnthalpy)),'SurfaceFlow:TotalEnthalpy','总焓与该条流线的入口不一致。');
    entropyScale=max([R,cp,abs(problem.referenceEntropy)]);
    assert(abs(entropy-problem.referenceEntropy)<=options.relativeStateTolerance*entropyScale, ...
        'SurfaceFlow:Entropy','熵与该条流线的入口不一致。');
end

%% 三、一般曲面坐标的二阶流线方程，不投影或裁剪积分状态
geometricSecond=B(:,1)*v(1)^2+2*B(:,2)*v(1)*v(2)+B(:,3)*v(2)^2;
connection=(C\(unitA'*geometricSecond))./lengths;
contravariantP=(C\(sample.pressureGradient./lengths))./lengths;
pressureTurn=(normT^2*contravariantP-v*(v'*sample.pressureGradient))/(rho*U^2);
rhs=[v;-connection-pressureTurn];
assert(all(isfinite(rhs)),'SurfaceFlow:Finite','流线导数非有限。');
info=struct('sample',sample,'position',sample.position,'tangent',tangent,'normal',normal, ...
    'coordinateUnitTangent',coordinateUnit,'unitTangentError',abs(normT-1), ...
    'totalEnthalpy',totalH,'entropy',entropy,'inputAlongMomentumResidualPaPerM',momentum, ...
    'inputAlongEnergyResidualJPerKgM',tangent'*energyGradient, ...
    'inputAlongEntropyResidualJPerKgKM',tangent'*entropyGradient, ...
    'thermalPropertyMode',thermal.thermalPropertyMode,'specificHeat',cp, ...
    'geodesicCurvatureVector',-(gradP-tangent*(tangent'*gradP))/(rho*U^2));
end
