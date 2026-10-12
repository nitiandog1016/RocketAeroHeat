function testContinuumRegime()
% 连续介质筛选：解析平均自由程、局部稀薄区及准流线下游传递。
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% 1. 标准状态平均自由程应为约 0.07 微米
lambda=airMeanFreePath(101325,300,1.8e-5,287);
expected=1.8e-5/101325*sqrt(pi*287*300/2);
assert(abs(lambda/expected-1)<1e-12);
assert(lambda>5e-8 && lambda<1e-7);
assert(isnan(airMeanFreePath(-1,300,1.8e-5,287)));

%% 2. 第二面元进入滑移区，其下游虽恢复局部压力仍不重新套用热流
m.F=repmat([1 2 3],4,1);
m.V=[0 0 0;1 0 0;0 1 0];
m.boundaryLayerEdge=struct('pressure',[101325;40;101325;101325], ...
    'temperature',300*ones(4,1), ...
    'dynamicViscosity',1.8e-5*ones(4,1), ...
    'algebraicallyValid',true(4,1),'method','解析测试外缘');
m.surfaceStreamlines=struct('predecessorFace',[0;1;2;3], ...
    'pathCost',(0:3)','connectedFace',true(4,1), ...
    'velocitySource','解析测试外缘');
fs=struct('pressure',101325,'temperature',300, ...
    'dynamicViscosity',1.8e-5,'R',287);
cfg.heating.continuumLength=0.01;
checked=assessContinuumRegime(m,fs,cfg);
r=checked.continuumRegime;
assert(r.regimeID(1)==1 && r.regimeID(2)==2);
assert(isequal(r.directUnresolvedFace,[false;true;false;false]));
assert(isequal(r.downstreamUnresolvedFace,[false;false;true;true]));
assert(isequal(r.applicabilityID,[1;2;5;5]));
assert(all(r.unresolvedFace(2:4)));
assert(contains(r.method,'平均自由程'));
tooLoose=cfg; tooLoose.heating.maxContinuumKnudsen=0.02;
failed=false;
try, assessContinuumRegime(m,fs,tooLoose); catch, failed=true; end
assert(failed);

%% 3. 自由来流本身稀薄时，即使局部压力升高也不能恢复连续介质
rare=fs; rare.pressure=40;
checked=assessContinuumRegime(m,rare,cfg);
assert(all(checked.continuumRegime.unresolvedFace));
assert(all(checked.continuumRegime.directUnresolvedFace));

%% 4. 默认模型高度上限附近，鼻端经验式须明确退出适用域
evalc('flightCfg=createCaseConfig();');
flightCfg.flight.altitude=80000;
high=calcFreestream(flightCfg.flight);
point.stagnation=struct('representativeFaceID',1, ...
    'representativePosition',[0 0 0], ...
    'hasPhysicalRadius',true,'physicalRadius',0.01, ...
    'radiusSource','测试给定物理半径','numberComponents',1);
point=calcStagnationHeating(point,high,flightCfg);
assert(~point.stagnationHeating.available);
assert(point.stagnationHeating.freestreamKnudsen>0.01);
assert(contains(point.stagnationHeating.status,'Kn'));
fprintf('Continuum regime passed: lambda_STP=%.3g m, Kn_80km=%.3g.\n', ...
    lambda,point.stagnationHeating.freestreamKnudsen);
end
