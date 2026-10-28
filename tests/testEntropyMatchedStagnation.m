function testEntropyMatchedStagnation()
%TESTENTROPYMATCHEDSTAGNATION 内部驻点工程质量/热力近似的独立局部校核。
% 模块作用：核对独立完全气体总压导数、耦合方程、有限s热力极限和
% 同剖面质量/热流尺度；本例没有制造质量厚度，也没有使用外部CFD。
% 不将鼻端局部通过当作鼻部/后体全域或实验热流验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
fs=makeState(220,5000,4);geo=struct('noseRadius',.1,'baseRadius',.5, ...
    'coneAngleDeg',25,'tipX',1);
cfg=struct('heating',struct('airThermalModel','caloricallyPerfect','Prandtl',.72));
cfg.externalFlow.entropyShock.stationCount=81;
result=solveEntropyMatchedStagnation(fs,geo,300,cfg);
assert(~result.available && result.diagnosticAvailable,result.status);
assert(isnan(result.heatFlux)&&isnan(result.massThickness));
estimate=result.conditionalEstimate;
assert(~result.defaultSolverEnabled && ~result.fullBodySolved && ~result.experimentallyValidated);
assert(~result.physicalBoundaryLayerClosed && ~result.entropyVorticityTermsSolved);
assert(result.diagnosticOnly && ~result.eligibleForSurfaceCoverage);
assert(strcmp(estimate.valueKind,'conditionalEstimate') && ~estimate.eligibleForSurfaceCoverage);
assert(result.localThermalMassClosureSolved && result.continuumScreeningPassed);
assert(result.stagnationMomentumRelativeMismatch>0 && ...
    result.stagnationMomentumRelativeMismatch>result.momentumApproximation.limit);
assert(abs(result.stagnationPressureParameter+.5- ...
    .5*result.stagnationMomentumRelativeMismatch)<1e-12);
assert(estimate.velocityGradient>estimate.fixedEntropyVelocityGradient);
assert(result.massClosureResidual<1e-8 && estimate.heatFlux>0);
assert(estimate.profile.massDefinitionResidual<1e-12);
assert(abs(estimate.massThickness/estimate.fixedEntropyProfile.massThickness- ...
    sqrt(estimate.fixedEntropyVelocityGradient/estimate.velocityGradient))<1e-8);
assert(abs(estimate.heatFlux/estimate.fixedEntropyProfile.heatFlux- ...
    sqrt(estimate.velocityGradient/estimate.fixedEntropyVelocityGradient))<1e-8);

%% 1. 独立数值完全气体公式导数，不调用被测热力内核作参考
w=1e-4;coarse=(p0Caloric(fs,w)-p0Caloric(fs,-w))/(2*w);
fine=(p0Caloric(fs,w/2)-p0Caloric(fs,-w/2))/w;
B=(4*fine-coarse)/3;
assert(abs(result.totalPressureAngleDerivative/B-1)<1e-8);
assert(abs(estimate.velocityGradient^2-estimate.fixedEntropyVelocityGradient^2- ...
    2*result.totalPressureMassDerivative*estimate.velocityGradient*estimate.massThickness)/ ...
    estimate.velocityGradient^2<1e-8);

%% 2. 有限小s独立热力恢复收敛到耦合应变率，不能保持固定驻点熵
ratio=[.01 .005 .0025];errors=zeros(size(ratio));
cp=fs.gamma*fs.R/(fs.gamma-1);T0=fs.temperature+fs.speed^2/(2*cp);
for j=1:numel(ratio)
    s=ratio(j)*geo.noseRadius;
    Psi=estimate.edgeDensity*estimate.velocityGradient*estimate.massThickness*s^2;
    w=4*result.shockTable.noseExpansion.A*Psi/(fs.density*fs.speed*geo.noseRadius^2);
    p0=p0Caloric(fs,w);
    pe=fs.pressure+(estimate.edgePressure-fs.pressure)*cos(s/geo.noseRadius)^2;
    T=T0*(pe/p0)^((fs.gamma-1)/fs.gamma);U=sqrt(2*cp*(T0-T));
    errors(j)=abs(U/(estimate.velocityGradient*s)-1);
end
assert(all(diff(errors)<0) && errors(end)<2e-5);
fprintf('Local engineering mass/profile coupling: strain %.9g vs fixed-entropy %.9g 1/s, q %.9g W/m2, mass %.9g m, finite-s relative error %.3g.\n', ...
    estimate.velocityGradient,estimate.fixedEntropyVelocityGradient,estimate.heatFlux,estimate.massThickness,errors(end));
fprintf('Leading-order outer momentum mismatch %.6g; source BL/vorticity closure remains unsolved.\n', ...
    result.stagnationMomentumRelativeMismatch);

%% 3. NASA7同热物性/同剖面耦合及域外/非零攻角拒绝
thermal=cfg;thermal.heating.airThermalModel='thermallyPerfect';
therm=solveEntropyMatchedStagnation(fs,geo,300,thermal);
assert(~therm.available && therm.diagnosticAvailable,therm.status);
assert(therm.massClosureResidual<1e-8 && therm.conditionalEstimate.heatFlux>0);
assert(therm.derivativeRelativeUncertainty<1e-4);
independentB=thermalDerivativeReference(fs);
assert(abs(therm.totalPressureAngleDerivative/independentB-1)<1e-7);
fprintf('NASA7 independent implicit normal-shock derivative error %.3g.\n', ...
    abs(therm.totalPressureAngleDerivative/independentB-1));
assert(therm.conditionalEstimate.profile.massDefinitionResidual<1e-12);
oblique=fs;oblique.alpha=1;refused=solveEntropyMatchedStagnation(oblique,geo,300,cfg);
assert(~refused.available && isnan(refused.heatFlux));
refused=solveEntropyMatchedStagnation(fs,geo,3000,thermal);
assert(~refused.available && isnan(refused.heatFlux));
strict=cfg;strict.massProfile.maximumMomentumMismatch=.005;
refused=solveEntropyMatchedStagnation(fs,geo,300,strict);
assert(~refused.available && isnan(refused.heatFlux));
assert(refused.stagnationMomentumRelativeMismatch>.005);
rarefied=makeState(220,.05,4);refused=solveEntropyMatchedStagnation(rarefied,geo,300,cfg);
assert(~refused.available && isnan(refused.heatFlux));
assert(~refused.continuumScreeningPassed && refused.continuumScreening.maximumKnudsen>=.01);
save(fullfile(root,'outputs','v52_validation','stagnation_mass_profile_coupling.mat'), ...
    'result','therm','errors','fs','geo');
fprintf('V52 STAGNATION MASS PROFILE COUPLING PASSED; internally computed mass thickness and explicit momentum mismatch, no physical BL/full-body/experimental claim.\n');
end
function fs=makeState(T,p,M)
R=287.05;g=1.4;fs=struct('temperature',T,'pressure',p,'R',R,'gamma',g,'Mach',M, ...
    'density',p/(R*T),'speed',M*sqrt(g*R*T),'alpha',0,'beta',0);
end
function p0=p0Caloric(fs,w)
g=fs.gamma;M=fs.Mach;mn2=M^2*(1-w);
r=(g+1)*mn2/((g-1)*mn2+2);
p2=fs.pressure*(1+2*g/(g+1)*(mn2-1));
T2=fs.temperature*(p2/fs.pressure)/r;
T0=fs.temperature*(1+(g-1)*M^2/2);p0=p2*(T0/T2)^(g/(g-1));
end
function B=thermalDerivativeReference(fs)
% 独立温度未知量正激波根 + 质量/动量/能量的隐式微分。
% 不复用被测热力核或其有限差分，T0由完整上游总焓反解且不随w变。
N=fs.speed;[~,h1]=airThermalProperties(fs.temperature);H0=h1+N^2/2;
T0=fzero(@(T)enthalpyResidual(T,H0),[fs.temperature 2500]);
T2=fzero(@(T)shockMomentum(T,fs,h1),[fs.temperature*1.0001,T0*(1-1e-6)]);
[cp2,h2,s2]=airThermalProperties(T2);[~,~,s0]=airThermalProperties(T0);
v2=sqrt(2*(H0-h2));rho2=fs.density*N/v2;p2=rho2*fs.R*T2;
p0=p2*exp((s0-s2)/fs.R);Nw=-N/2;
coeff=[cp2 v2;N*fs.R/v2 N*(1-fs.R*T2/v2^2)];
derivative=coeff\[N*Nw;-(fs.R*T2/v2+v2-2*N)*Nw];
Tw=derivative(1);vw=derivative(2);
B=p0*(Nw/N+(1-cp2/fs.R)*Tw/T2-vw/v2);
end
function value=shockMomentum(T,fs,h1)
[~,h]=airThermalProperties(T);v=sqrt(fs.speed^2+2*(h1-h));
rho=fs.density*fs.speed/v;
value=rho*fs.R*T+rho*v^2-fs.pressure-fs.density*fs.speed^2;
end
function value=enthalpyResidual(T,target)
[~,h]=airThermalProperties(T);value=h-target;
end
