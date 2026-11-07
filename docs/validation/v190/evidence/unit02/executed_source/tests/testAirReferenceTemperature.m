function report=testAirReferenceTemperature()
% 独立积分cp核对参考焓、常cp极限和反馈导数，不只复写实现公式。
fs=struct('R',287.05,'gamma',1.4); cfg=struct;
Te=[180;350;1200]; Tw=[300;800;1800]; Taw=[600;1600;2200];
[Tr,state]=airReferenceTemperature(Te,Tw,Taw,fs,cfg);
residual=zeros(3,1);
for j=1:3
    % 分段积分绕开NASA系数1000K处舍入不连续；该差异单列。
    lhs=integrateCp(Te(j),Tr(j));
    rhs=0.5*integrateCp(Te(j),Tw(j))+0.22*integrateCp(Te(j),Taw(j));
    residual(j)=abs(lhs-rhs);
end
fprintf('Independent cp integral residuals J/kg: %.17g %.17g %.17g\n',residual);
% 此门槛事前允许原NASA分段圆整差，不放宽参考焓反解自身的严格门槛。
assert(max(residual)<0.01);
oldCfg.heating.referencePropertyModel='linearTemperature';
linear=airReferenceTemperature(Te,Tw,Taw,fs,oldCfg);
assert(max(abs(Tr-linear))>10);
constant.heating.airThermalModel='caloricallyPerfect';
exact=airReferenceTemperature(Te,Tw,Taw,fs,constant);
assert(max(abs(exact-linear))<1e-10);
step=0.01;
plus=airReferenceTemperature(Te,Tw+step,Taw,fs,cfg);
minus=airReferenceTemperature(Te,Tw-step,Taw,fs,cfg);
derivativeError=max(abs((plus-minus)/(2*step)-state.temperatureDerivativeByWallTemperature));
assert(derivativeError<1e-7);
bad=airReferenceTemperature([70;300;NaN],[300;2600;300],[400;800;400],fs,cfg);
assert(all(isnan(bad)));
failed=false;
try, airReferenceTemperature(300,[300;400],800,fs,cfg); catch, failed=true; end
assert(failed);
failed=false; invalid=cfg; invalid.heating.referencePropertyModel='made-up';
try, airReferenceTemperature(300,300,800,fs,invalid); catch, failed=true; end
assert(failed);
% 相同基准壁温下h应保持，升温导数与中心差分一致。
p=2500*ones(3,1); mu=1.458e-6*Tr.^1.5./(Tr+110.4);
law=struct('edgeTemperature',Te,'edgePressure',p,'exponent',[.5;.2;.5], ...
    'referenceViscosity0',mu.*fs.R.*Tr./p,'gasConstant',fs.R, ...
    'thermalPropertyMode','thermallyPerfect','referencePropertyModel','referenceEnthalpy', ...
    'referenceCp0',airEnthalpySecantCp(Taw,Tw));
hBase=[100;200;300]; [h,dh]=evaluateWallConvection(hBase,Taw,Tw,law);
assert(max(abs(h-hBase))<1e-9);
hp=evaluateWallConvection(hBase,Taw,Tw+step,law);
hm=evaluateWallConvection(hBase,Taw,Tw-step,law);
feedbackDerivativeError=max(abs((hp-hm)/(2*step)-dh));
assert(feedbackDerivativeError<1e-6);
% 保留旧参考模式，历史law没有标记时仍按线性温度反馈。
legacy=rmfield(law,'referencePropertyModel');
oldMu=1.458e-6*linear.^1.5./(linear+110.4);
legacy.referenceViscosity0=oldMu.*fs.R.*linear./p;
assert(max(abs(evaluateWallConvection(hBase,Taw,Tw,legacy)-hBase))<1e-9);
testBoundaryLayerState(); testAirThermalProperties(); testWallThermalResponse();
report=struct('passed',true,'cpIntegralResidualJkg',residual, ...
    'maximumChangeFromLinearK',max(abs(Tr-linear)), ...
    'temperatureDerivativeError',derivativeError, ...
    'feedbackDerivativeError',feedbackDerivativeError, ...
    'newPhysicalQualification',false);
disp(report);
end

function value=integrateCp(a,b)
% AbsTol不足以约束百万量级积分；显式给出相对精度与原分段点。
points=[150 1000]; points=points(points>min(a,b)&points<max(a,b));
value=integral(@(T)airThermalProperties(T),a,b, ...
    'AbsTol',1e-7,'RelTol',1e-11,'Waypoints',points);
end
