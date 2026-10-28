function report=testNonsimilarLaminarModule(folder)
%TESTNONSIMILARLAMINARMODULE 正式条件模块的可移植定向数值与合同测试。
% 模块一：闭式Couette剪切/黏性发热，核对热流与剪切符号和无源方程。
% 模块二：明确制造源的非零流向场，核对不均匀n、不等步BDF1/2和Newton。
% 模块三：NASA7焓/导热身份、无效剖面与公共适用声明拒绝、失败保留NaN。
% 模块四：保存最终源码SHA256与真实MAT/JSON；不运行长球鼻study，
% 不增加生产覆盖，不用制造场或自报适用声明冒充实验/物理解。
root=fileparts(fileparts(mfilename('fullpath')));
library=fullfile(root,'regionalFlow','nonsimilar');
addpath(root,library);
if nargin<1,folder=fullfile(root,'outputs','v61_validation','formal_module');end
if ~exist(folder,'dir'),mkdir(folder);end
gas=struct('R',287.05287,'gamma',1.4,'thermalModel','caloricallyPerfect', ...
    'Prandtl',.72,'transport','constant','constantViscosity',1.8e-5);
cp=gas.gamma*gas.R/(gas.gamma-1);N=41;n=linspace(0,.002,N)';

%% 模块一：无源Couette闭式解，检查剪切、黏性功及壁热流
A=8e4;B=2e4;Tw=300;
c=station(.05,.1,n,A*n,Tw+B*n-gas.Prandtl*A^2*n.^2/(2*cp),gas);
h=c;h.s=.049;
residual=nonsimilarLaminarResidual(c,h,gas);
k=gas.constantViscosity*cp/gas.Prandtl;
couette=struct('maxMomentumResidual',max(abs(residual.momentum)), ...
    'maxEnergyResidual',max(abs(residual.energy)), ...
    'wallHeatFlux',residual.wallHeatFluxIntoWall,'analyticWallHeatFlux',k*B, ...
    'wallShear',residual.wallShearStress,'analyticWallShear',gas.constantViscosity*A, ...
    'continuity',max(abs(residual.continuityIntervalBalance)));
assert(couette.maxMomentumResidual<1e-6 && couette.maxEnergyResidual<2e-5);
assert(abs(couette.wallHeatFlux-k*B)<1e-7);
assert(abs(couette.wallShear-gas.constantViscosity*A)<1e-10);
assert(~residual.forcingPresent && ~residual.physicalSolution && ~residual.eligibleForSurfaceCoverage);
domain=struct('axisymmetric',true,'angleOfAttackDeg',0,'sideslipDeg',0, ...
    'flowRegime','continuum','layerState','laminar','crossflowPresent',false, ...
    'separationPresent',false,'thinLayerVerified',true, ...
    'inletIdentity','闭式Couette验证剖面；不是球鼻物理入口', ...
    'wallTemperatureHistoryIdentity','闭式Couette全段均匀Tw=300K给定');
wrapped=solveNonsimilarLaminarStation(c,h,gas,domain);
assert(wrapped.available && wrapped.numericalConverged);
assert(~wrapped.physicalSolution && ~wrapped.eligibleForSurfaceCoverage);
assert(~wrapped.domainDeclarationAutomaticallyVerified && isequal(wrapped.domainDeclaration,domain));
assert(isnan(wrapped.wallTemperatureSolved) && wrapped.wallTemperatureInput==Tw);

%% 模块二：不均匀物理n和不等s步长的非零流向制造场
% u=(A0+A1*s)n，T固定，v=-A1*n^2/2；给定显式动量/能量源。
% 平面定r可精确再现BDF2；这项测试仅为迁移回归，不追加制造场研究。
n=n(end)*(linspace(0,1,N)').^1.4;
A0=8e4;A1=2e5;rho=5000/(gas.R*300);
c=manufactured(.05);h1=manufactured(.049);h2=manufactured(.0475);
a=A0+A1*c.s;
fm=.5*rho*a*A1*n.^2;
fe=.5*rho*a^2*A1*n.^3-gas.constantViscosity*a^2;
forcing=struct('momentum',fm,'energy',fe);
rr=nonsimilarLaminarResidual(c,[h1,h2],gas,forcing);
firstOrder=nonsimilarLaminarResidual(c,h1,gas,forcing);
ds=c.s-h1.s;oldDs=h1.s-h2.s;
expectedWeights=[(2*ds+oldDs)/(ds*(ds+oldDs)), ...
    -(ds+oldDs)/(ds*oldDs),ds/(oldDs*(ds+oldDs))];
exact=struct('maxMomentumResidual',max(abs(rr.momentum)), ...
    'maxEnergyResidual',max(abs(rr.energy)), ...
    'maximumNormalVelocityError',max(abs(rr.normalVelocity+.5*A1*n.^2)), ...
    'streamwiseMomentumMagnitude',max(abs(rr.momentumConvection)), ...
    'streamwiseEnergyMagnitude',max(abs(rr.energyConvection)), ...
    'normalVelocityAtTop',rr.normalVelocity(end), ...
    'nonuniformNormalSpacingRatio',max(diff(n))/min(diff(n)), ...
    'streamwiseHistoryStepSizesM',[ds,oldDs]);
assert(exact.maxMomentumResidual<1e-6 && exact.maxEnergyResidual<3e-4);
assert(exact.maximumNormalVelocityError<1e-10);
assert(exact.streamwiseMomentumMagnitude>1 && exact.streamwiseEnergyMagnitude>1);
assert(abs(exact.normalVelocityAtTop)>.01 && exact.nonuniformNormalSpacingRatio>2);
assert(rr.streamwiseDifferenceOrder==2 && firstOrder.streamwiseDifferenceOrder==1);
assert(max(abs(rr.streamwiseWeights-expectedWeights))<1e-10);
assert(max(abs(firstOrder.streamwiseWeights-[1,-1]/ds))<1e-10);
assert(max(abs(rr.velocityStreamwiseDerivative-A1*n))<1e-7);
assert(max(abs(firstOrder.velocityStreamwiseDerivative-A1*n))<1e-7);
assert(~rr.normalTopVelocityPrescribed && rr.forcingPresent);
guess=c;guess.u=c.u.*(1+.03*sin(pi*n/n(end)));
guess.T=c.T+3*sin(pi*n/n(end));
advanced=advanceNonsimilarLaminarStation(guess,[h1,h2],gas,struct,forcing);
assert(advanced.available && advanced.converged && advanced.residualHistory(end)<1e-9);
exact.newtonIterations=advanced.iterations;
exact.newtonResidualHistory=advanced.residualHistory;
exact.maxAdvancedVelocityError=max(abs(advanced.profile.u-c.u));
exact.maxAdvancedTemperatureError=max(abs(advanced.profile.T-c.T));
assert(exact.maxAdvancedVelocityError<2e-5 && exact.maxAdvancedTemperatureError<2e-5);
assert(~advanced.physicalSolution && ~advanced.eligibleForSurfaceCoverage);

%% 模块三：NASA7与拒绝合同；无法求解时保留未知
ng=gas;ng.thermalModel='thermallyPerfect';ng.transport='sutherland';
c=station(.05,.1,n,A*n,300+B*n,ng);h=c;h.s=.049;
rr=nonsimilarLaminarResidual(c,h,ng);[cc,hh]=airThermalProperties(c.T);
nasa=struct('enthalpyMaximumDifference',max(abs(rr.enthalpy-hh)), ...
    'heatCapacityMaximumDifference',max(abs(rr.heatCapacity-cc)), ...
    'totalEnthalpyIdentityError',max(abs(rr.totalEnthalpy-hh-.5*c.u.^2)), ...
    'conductivityIdentityError',max(abs(rr.conductivity-rr.viscosity.*cc/ng.Prandtl)));
assert(nasa.enthalpyMaximumDifference==0 && nasa.heatCapacityMaximumDifference==0);
assert(nasa.totalEnthalpyIdentityError<1e-9 && nasa.conductivityIdentityError==0);
rejections=struct;
bad=h;bad.n(end)=bad.n(end)+1e-8;
rejections.differentPhysicalGrid=refused(@()nonsimilarLaminarResidual(c,bad,ng));
bad=c;bad.s=0;
rejections.singularStagnationStation=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=c;bad.T(3)=79;
rejections.outsideNASAProperty=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=c;bad.He=bad.He+100;
rejections.inconsistentEdgeEnthalpy=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=c;bad.sources=rmfield(bad.sources,'profile');
rejections.missingProfileIdentity=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=c;bad.u(3)=-1;
rejections.reverseCurrentVelocity=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=h;bad.u(3)=-1;
rejections.reverseHistoryVelocity=refused(@()nonsimilarLaminarResidual(c,bad,ng));
bad=c;bad.u=bad.u+1i;
rejections.complexCurrentVelocity=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=c;bad.T=bad.T+1i;
rejections.complexCurrentTemperature=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=c;bad.n=bad.n+1i;
rejections.complexPhysicalGrid=refused(@()nonsimilarLaminarResidual(bad,h,ng));
bad=ng;bad.gamma=Inf;
rejections.infiniteGasGamma=refused(@()nonsimilarLaminarResidual(c,h,bad));
badf=struct('energy',ones(N,1)*1i);
rejections.complexManufacturedSource=refused(@()nonsimilarLaminarResidual(c,h,ng,badf));
bad=h;bad.s=c.s-eps(c.s);bad.r=1e300;
rejections.finiteInputOverflow=refused(@()nonsimilarLaminarResidual(c,bad,ng));
rejections.unknownNewtonOption=refused(@()advanceNonsimilarLaminarStation(c,h,ng,struct('inventedTolerance',1)));
stopped=advanceNonsimilarLaminarStation(c,h,ng,struct('maximumIterations',0));
rejections.workLimitLeavesUnknown=~stopped.available && isnan(stopped.wallHeatFluxIntoWall) && isempty(stopped.profile);
% 这些是物性边界可执行性检查，不宣称一次迭代已收敛。
hot=station(.05,.1,n,A*n,2500+zeros(N,1),gas);oldHot=hot;oldHot.s=.049;
hot.u=hot.u.*(1+.001*sin(pi*n/n(end)));
hotResult=advanceNonsimilarLaminarStation(hot,oldHot,gas,struct('maximumIterations',1));
rejections.upperPropertyBoundJacobianHasOneSidedDifference=isstruct(hotResult);
hot.T(2:end-1)=1500;oldHot.T=hot.T;
hotResult=advanceNonsimilarLaminarStation(hot,oldHot,gas, ...
    struct('maximumIterations',1,'finiteDifferenceStep',.9));
rejections.largeJacobianDifferenceStaysInsideBothBounds=isstruct(hotResult);
bad=domain;bad.angleOfAttackDeg=20;
rejections.obliqueDomain=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.crossflowPresent=true;
rejections.crossflowDomain=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.separationPresent=true;
rejections.separatedDomain=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.thinLayerVerified=false;
rejections.uncheckedThinLayer=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.wallTemperatureHistoryIdentity='';
rejections.missingWallHistory=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.inletIdentity=' ';
rejections.missingInletIdentity=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.flowRegime='rarefied';
rejections.rarefiedDomain=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.layerState='turbulent';
rejections.turbulentDomain=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.axisymmetric=false;
rejections.nonaxisymmetricDomain=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.sideslipDeg=1;
rejections.sideslipDomain=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.ignoreEnergy=true;
rejections.unknownDomainField=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
bad=domain;bad.thinLayerVerified=1;
rejections.numericTruthNotLogicalDeclaration=refused(@()solveNonsimilarLaminarStation(c,h,ng,bad));
assert(all(structfun(@(q)logical(q),rejections)));

%% 模块四：最终源码绑定，保留实际剖面与终态证据
files={'regionalFlow/nonsimilar/nonsimilarLaminarResidual.m', ...
    'regionalFlow/nonsimilar/advanceNonsimilarLaminarStation.m', ...
    'regionalFlow/nonsimilar/solveNonsimilarLaminarStation.m', ...
    'tests/testNonsimilarLaminarModule.m','airThermalProperties.m'};
bindings=repmat(struct('path','','sha256',''),numel(files),1);
for j=1:numel(files)
    bindings(j).path=files{j};bindings(j).sha256=sha256(fullfile(root,files{j}));
end
report=struct('dateUTC',char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss'Z'")), ...
    'tested',true,'testPassed',true,'physicalSolution',false, ...
    'eligibleForSurfaceCoverage',false,'newSolvedAreaM2',0, ...
    'scope','正式条件模块的单站定向数值/合同回归；非球鼻实验或新增生产物理解', ...
    'couette',couette,'nonzeroStreamwiseManufactured',exact, ...
    'NASA7Contract',nasa,'rejections',rejections,'rejectionCheckCount',numel(fieldnames(rejections)), ...
    'sourceBindings',bindings);
save(fullfile(folder,'nonsimilar_module_evidence.mat'),'report','wrapped','advanced','c','h','rr','firstOrder','stopped');
fid=fopen(fullfile(folder,'nonsimilar_module_summary.json'),'w','n','UTF-8');
assert(fid>=0);clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear clean;
fprintf('PASS: nonsimilar module numeric/contracts, %d rejection checks, source-bound evidence.\n',report.rejectionCheckCount);

    function z=manufactured(s)
        z=station(s,.1,n,(A0+A1*s)*n,300+zeros(N,1),gas);
    end
end

function z=station(s,r,n,u,T,g)
% 辅助模块：明确剖面来源，生成同物性的外缘总焓身份。
if strcmp(g.thermalModel,'thermallyPerfect'),[~,he]=airThermalProperties(T(end));
else,he=g.gamma*g.R/(g.gamma-1)*T(end);end
z=struct('s',s,'r',r,'p',5000,'p_s',0,'n',n,'u',u,'T',T, ...
    'Ue',u(end),'Te',T(end),'He',he+.5*u(end)^2,'Tw',T(1), ...
    'sources',struct('geometry','明确合成测试几何', ...
    'pressure','闭式验证场给定压力，不是球鼻外流', ...
    'profile','闭式连续u(s,n)/T(s,n)点值剖面', ...
    'edge','同验证场/同物性外界','wallTemperature','给定验证壁温'));
end

function yes=refused(call)
% 辅助模块：只把预期的合同异常算拒绝；其他异常重新抛出。
yes=false;try,call();catch exception
    if startsWith(exception.identifier,'RocketAeroHeat:NonSimilar:'),yes=true;
    else,rethrow(exception);end
end
end

function digest=sha256(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
