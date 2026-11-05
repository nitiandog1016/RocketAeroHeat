function freestream=calcSU2Freestream(cfg)
%CALCSU2FREESTREAM 内嵌计算页按显式物性模型建立来流参考。
% 模块一。沿用坐标、静温/静压和输入Mach，选择已验证的理想空气核。
% 模块二。同步更新声速、速度、密度和动压，避免UI使用旧定比热速度。
% 模块三。总状态由同一焓/熵核反解，定义为不含SST k的热力学停滞状态。
freestream=calcFreestream(cfg.flight);
if ~isfield(cfg.su2,'thermodynamicModel') || ~strcmp(cfg.su2.thermodynamicModel,'idealAir'),return;end
air=idealAirThermodynamicState('PT',freestream.pressure,freestream.temperature);
speed=freestream.Mach*air.sound;
total=idealAirThermodynamicState('hs',air.h+speed^2/2,air.s);
muRef=1.458e-6*273.15^1.5/(273.15+110.4);
mu=muRef*(air.T/273.15)^1.5*383.55/(air.T+110.4);
freestream.density=air.rho;freestream.speedOfSound=air.sound;
freestream.R=air.gasConstant;freestream.gamma=air.gamma;freestream.dynamicViscosity=mu;
freestream.speed=speed;freestream.Vinf=speed;
freestream.flowVelocityVector=speed*freestream.unitFlowDirection;
freestream.flightVelocityVector=speed*freestream.unitFlightDirection;
freestream.dynamicPressure=.5*air.rho*speed^2;freestream.qinf=freestream.dynamicPressure;
freestream.reynoldsPerMeter=air.rho*speed/mu;freestream.kineticEnergyPerMass=speed^2/2;
freestream.totalTemperature=total.T;freestream.totalPressure=total.p;
freestream.thermodynamicModel=air.model;freestream.energyEntropyReference=air.reference;
freestream.thermalEnthalpy=air.h;freestream.thermalInternalEnergy=air.e;
freestream.referenceConvention='preserve-Mach-static-PT-direction';
freestream.atmosphere.density=air.rho;freestream.atmosphere.speedOfSound=air.sound;
freestream.atmosphere.R=air.gasConstant;freestream.atmosphere.gamma=air.gamma;
freestream.atmosphere.dynamicViscosity=mu;
end
