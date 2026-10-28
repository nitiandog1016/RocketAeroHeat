function testAxisymmetricLaminarMassProfile()
%TESTAXISYMMETRICLAMINARMASSPROFILE 局部相似剖面的独立极限及拒绝校核。
% 模块作用：验证Blasius/Falkner-Skan/轴对称驻点的独立常密度基准，
% Pr=1可压缩平板的Crocco焓关系，同剖面厚度、热流和真实域外状态。
% 不据此声称真实球锥耦合、生产热流或公开实验数据已经验证。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
fs=struct('R',287.05,'gamma',1.4);
caloric=struct('heating',struct('airThermalModel','caloricallyPerfect','Prandtl',1));
cp=fs.gamma*fs.R/(fs.gamma-1);

%% 1. 低Mach常密度极限：独立Blasius及Falkner-Skan壁面剪切基准
% 基准来自SciPy的独立三状态ODE，不调用本模块建立参考。
benchBeta=[0 .2 .5]; benchSlope=[.46959998836103 .686708181032189 .927680039836679];
for j=1:numel(benchBeta)
    [e,s]=makeStation(300,101325,1e-3,benchBeta(j),fs,cp);
    r=solveAxisymmetricLaminarMassProfile(fs,e,300,s,caloric);
    assert(r.available,r.status);
    assert(abs(r.wallDimensionlessVelocitySlope/benchSlope(j)-1)<2e-6);
    assert(r.localSimilarity && ~r.nonSimilarTermsSolved && ~r.massMatchingSolved);
    assert(r.etaEdge<r.numericalOptions.etaMax); % 实际外缘必须在数值边界以内
    assert(r.inputMomentumRelativeResidual<1e-12);
    assert(r.massDefinitionResidual<1e-12);
    assert(r.delta>r.massThickness && r.massThickness>r.thetaPhysical);
    assert(abs(r.profile.temperature(end)-300)<1e-8);
    if j==1
        thetaInfinite=.46959998836103;
        assert(abs(r.dimensionlessIntegrals.momentum/thetaInfinite-1)<.003);
        % 有限边缘的质量积分独立由f(eta_edge)给出，不能以5.27代替。
        ref=independentSimilarity(benchBeta(j),r.etaEdge,benchSlope(j));
        assert(abs(r.dimensionlessIntegrals.mass/ref(1)-1)<2e-6);
        assert(abs(r.dimensionlessIntegrals.momentum/ref(4)-1)<2e-6);
        assert(abs(r.massThickness/r.thetaPhysical-5.27)>.05);
    end
    fprintf('Constant-density beta %.2f: fpp0 %.12g, independent %.12g, edge %.6g.\n', ...
        benchBeta(j),r.wallDimensionlessVelocitySlope,benchSlope(j),r.etaEdge);
end

%% 2. Pr=1可压缩平板：Crocco关系及物理坐标积分、壁面导数
[e,s]=makeStation(400,8000,1000,0,fs,cp);
r=solveAxisymmetricLaminarMassProfile(fs,e,300,s,caloric);
assert(r.available,r.status);
F=r.profile.states(2,:); g=r.profile.states(4,:);
gw=cp*300/e.totalEnthalpy;
assert(max(abs(g-(gw+(1-gw)*F)))<2e-6);
expectedStatic=(e.totalEnthalpy*g-e.speed^2*F.^2/2)/cp;
assert(max(abs(r.profile.temperature-expectedStatic))<1e-9);
assert(max(abs(r.profile.density.*fs.R.*r.profile.temperature/e.pressure-1))<1e-12);
% 用物理坐标上的rho*u积分进行独立复算，不复用质量厚度变换公式。
y=r.profile.physicalY;
massIndependent=trapz(y,(r.profile.density/e.density).*F);
thetaIndependent=trapz(y,(r.profile.density/e.density).*F.*(1-F));
deltaStarIndependent=trapz(y,1-(r.profile.density/e.density).*F);
assert(abs(massIndependent/r.massThickness-1)<2e-5);
assert(abs(thetaIndependent/r.thetaPhysical-1)<2e-5);
assert(abs(deltaStarIndependent-r.deltaStar)/r.delta<2e-5);
muWall=1.458e-6*300^1.5/(300+110.4); rhoWall=e.pressure/(fs.R*300);
expectedQ=muWall*rhoWall/e.density/r.thicknessScale* ...
    e.totalEnthalpy*r.wallDimensionlessEnthalpySlope;
expectedTau=muWall*rhoWall/e.density/r.thicknessScale* ...
    e.speed*r.wallDimensionlessVelocitySlope;
assert(abs(r.heatFlux/expectedQ-1)<1e-10 && r.heatFlux>0);
assert(r.heatFluxNumericallyResolved && r.fieldMetadata.heatFlux.numericallyResolved);
assert(abs(r.wallShearStress/expectedTau-1)<1e-10);
assert(strcmp(r.fieldMetadata.wallTemperatureInput.valueKind,'prescribed'));
assert(strcmp(r.fieldMetadata.temperature.valueKind,'input'));
assert(strcmp(r.fieldMetadata.massThickness.valueKind,'computed') && r.fieldMetadata.massThickness.valid);
assert(r.edgeChecks.velocityDefect<=.005+1e-10);
assert(r.edgeChecks.totalEnthalpyDefect<=.005+1e-10);
assert(r.edgeChecks.relativeShearGradient<=.005+1e-10);
assert(r.edgeChecks.relativeHeatGradient<=.005+1e-10);
assert(r.edgeChecks.heatGradientCriterionDefect<=1+1e-8);
assert(r.edgeChecks.relativeTemperatureDefect<=.005+1e-10);
assert(r.edgeChecks.relativeDensityDefect<=.005+1e-10);
fprintf('Compressible Pr1: q %.9g W/m2, tau %.9g Pa, mass %.9g m, theta %.9g m.\n', ...
    r.heatFlux,r.wallShearStress,r.massThickness,r.thetaPhysical);

%% 3. NASA7焓零点、恢复关系及解析驻点极限
thermal=struct('heating',struct('Prandtl',1));
[~,he]=airThermalProperties(e.temperature); e.totalEnthalpy=he+e.speed^2/2;
nasa=solveAxisymmetricLaminarMassProfile(fs,e,300,s,thermal);
assert(nasa.available,nasa.status);
[~,hw]=airThermalProperties(300);
physicalTotal=nasa.profile.totalEnthalpy;
expected=hw+(e.totalEnthalpy-hw)*nasa.profile.states(2,:);
assert(max(abs(physicalTotal-expected))/abs(e.totalEnthalpy-hw)<2e-6);
[~,hProfile]=airThermalProperties(nasa.profile.temperature);
assert(max(abs(physicalTotal-hProfile-nasa.profile.speed.^2/2))<1e-5);
assert(nasa.enthalpyOffset~=0);
% 在U=0低温下原NASA7总焓为负，移位只处理数值变量，不伪造物理焓。
[es,ss]=makeStagnation(100,1000,50,fs,cp,'thermallyPerfect');
assert(es.totalEnthalpy<0);
st=solveAxisymmetricLaminarMassProfile(fs,es,100,ss,thermal);
assert(st.available,st.status);
mu=1.458e-6*100^1.5/(100+110.4);
assert(abs(st.thicknessScale/sqrt(mu/es.density/(2*50))-1)<1e-12);
assert(st.wallShearStress==0 && st.massFluxPerMetric==0);
assert(abs(st.heatFlux)<1e-5 && st.massThickness>0);
assert(st.stagnationLimit && abs(st.wallDimensionlessVelocitySlope/benchSlope(3)-1)<2e-6);
assert(~st.stagnationPressureCurvatureVerified);
assert(abs(st.assumedStagnationPressureCurvature/(-es.density*50^2)-1)<1e-12);
% 有限非零站向同一轴对称驻点局部展开收敛。
[es2,ss2]=makeStation(300,101325,.01,.5,fs,cp);
es2.totalEnthalpy=cp*300+.01^2/2;
finite=solveAxisymmetricLaminarMassProfile(fs,es2,300,ss2,caloric);
[es0,ss0]=makeStagnation(300,101325,ss2.velocityGradient,fs,cp,'caloricallyPerfect');
limit=solveAxisymmetricLaminarMassProfile(fs,es0,300,ss0,caloric);
assert(finite.available && limit.available);
assert(abs(finite.massThickness/limit.massThickness-1)<1e-6);
assert(finite.edgeChecks.heatGradientMaximum<finite.numericalOptions.edgeAbsoluteGradientTolerance);
assert(finite.edgeChecks.heatGradientCriterionDefect<=1+1e-8);
assert(~finite.heatFluxNumericallyResolved); % 不把低于数值分辨率的热流声称有相对精度
fprintf('NASA7/negative-enthalpy stagnation checks passed.\n');

%% 4. 高动能有限外缘必须恢复静温/密度，不能只看总焓相对缺陷
[eh,sh]=makeStation(400,8000,1800,0,fs,cp);
[~,heh]=airThermalProperties(eh.temperature); eh.totalEnthalpy=heh+eh.speed^2/2;
highEnergy=struct('heating',struct('Prandtl',.72));
hr=solveAxisymmetricLaminarMassProfile(fs,eh,300,sh,highEnergy);
assert(hr.available,hr.status);
assert(hr.edgeChecks.relativeTemperatureDefect<=hr.numericalOptions.edgeTemperatureTolerance+1e-10);
assert(hr.edgeChecks.relativeDensityDefect<=hr.numericalOptions.edgeDensityTolerance+1e-10);
assert(abs(hr.profile.temperature(end)/eh.temperature-1)<=.005+1e-10);
assert(abs(hr.profile.density(end)/eh.density-1)<=.005+1e-10);
fprintf('High-energy edge: eta %.6g, T defect %.6g, rho defect %.6g.\n', ...
    hr.etaEdge,hr.edgeChecks.relativeTemperatureDefect,hr.edgeChecks.relativeDensityDefect);
% 扩大纯数值截断、保持相同物理边缘规则；边缘及壁面场量应收敛。
extended=highEnergy; extended.massProfile=struct('etaMax',16);
hrExtended=solveAxisymmetricLaminarMassProfile(fs,eh,300,sh,extended);
assert(hrExtended.available,hrExtended.status);
assert(abs(hrExtended.massThickness/hr.massThickness-1)<1e-5);
assert(abs(hrExtended.heatFlux/hr.heatFlux-1)<1e-5);
assert(abs(hrExtended.wallShearStress/hr.wallShearStress-1)<1e-5);
fprintf('etaMax12/16 convergence: mass %.3g, q %.3g, tau %.3g relative differences.\n', ...
    hrExtended.massThickness/hr.massThickness-1,hrExtended.heatFlux/hr.heatFlux-1, ...
    hrExtended.wallShearStress/hr.wallShearStress-1);

%% 5. 拒绝不兼容梯度、温度域外及未解析有限外缘，不能补零
bad=s; bad.pressureGradient=0; bad.velocityGradient=100;
bad.beta=2*bad.xi*bad.velocityGradient/(e.density*sutherland(e.temperature)* ...
    e.speed^2*bad.streamlineMetric^2);
% beta在本模块范围内；pe_s=0但U_s非零时动量残差严格拒绝。
rejected=solveAxisymmetricLaminarMassProfile(fs,e,300,bad,thermal);
assert(~rejected.available && rejected.statusCode==2);
assert(rejected.inputMomentumRelativeResidual==1 && rejected.inputMomentumResidual>0);
assert(rejected.inputMomentumDenominator>0);
checkBlank(rejected);
bad=s; bad.beta=.2;
rejected=solveAxisymmetricLaminarMassProfile(fs,e,300,bad,thermal);
assert(~rejected.available && rejected.statusCode==7); checkBlank(rejected);
rejected=solveAxisymmetricLaminarMassProfile(fs,e,70,s,thermal);
assert(~rejected.available && rejected.statusCode==1); checkBlank(rejected);
short=thermal; short.massProfile=struct('etaMax',4,'edgeBuffer',1);
rejected=solveAxisymmetricLaminarMassProfile(fs,e,300,s,short);
assert(~rejected.available && rejected.statusCode==5); checkBlank(rejected);
hot=e; [~,hhot]=airThermalProperties(2500);
hot.temperature=2500; hot.density=hot.pressure/(fs.R*2500);
hot.totalEnthalpy=hhot+hot.speed^2/2;
rejected=solveAxisymmetricLaminarMassProfile(fs,hot,2500,s,thermal);
assert(~rejected.available); checkBlank(rejected);
missing=e; missing.sources=rmfield(missing.sources,'pressure');
didError=false;
try, solveAxisymmetricLaminarMassProfile(fs,missing,300,s,thermal);
catch ex, didError=strcmp(ex.identifier,'RocketAeroHeat:LaminarMassProfile:Source'); end
assert(didError);
fprintf('Axisymmetric local-similarity mass profile tests passed; no geometry/experimental claim.\n');
end

function [e,s]=makeStation(T,p,U,beta,fs,cp)
e=struct('temperature',T,'pressure',p,'density',p/(fs.R*T), ...
    'speed',U,'totalEnthalpy',cp*T+U^2/2,'sources',struct);
for key={'temperature','pressure','density','speed','totalEnthalpy'}
    e.sources.(key{1})='独立剖面极限校核输入，不是几何外流验证';
end
metric=.1; mu=sutherland(T); xi=e.density*mu*U*metric^2*.5;
dU=beta*e.density*mu*U^2*metric^2/(2*xi);
s=struct('streamlineMetric',metric,'xi',xi,'beta',beta, ...
    'pressureGradient',-e.density*U*dU,'velocityGradient',dU,'sources',struct);
for key={'streamlineMetric','xi','beta','pressureGradient','velocityGradient','wallTemperature'}
    s.sources.(key{1})='独立解析/数值输入；未用于真实球锥';
end
end
function [e,s]=makeStagnation(T,p,a,fs,cp,mode)
[e,s]=makeStation(T,p,1,.5,fs,cp);
e.speed=0;
if strcmp(mode,'thermallyPerfect'), [~,e.totalEnthalpy]=airThermalProperties(T);
else, e.totalEnthalpy=cp*T; end
s.streamlineMetric=0; s.xi=0; s.velocityGradient=a; s.pressureGradient=0;
s.stagnationVelocityGradient=a;
s.sources.stagnationVelocityGradient='轴对称驻点解析输入a；未由球锥耦合求得';
end
function y=independentSimilarity(beta,edge,slope)
% 独立三状态初值ODE，再同时积分物理常密度动量积分，不复用BVP。
opt=odeset('RelTol',2e-10,'AbsTol',1e-12);
[~,states]=ode45(@(~,z)[z(2);z(3);-z(1)*z(3)+beta*(z(2)^2-1); ...
    z(2)*(1-z(2))],[0 edge],[0;0;slope;0],opt);
y=states(end,:);
end
function checkBlank(r)
assert(isnan(r.massThickness) && isnan(r.thetaPhysical) && isnan(r.heatFlux));
assert(~r.fieldMetadata.massThickness.valid && ~r.fieldMetadata.heatFlux.valid);
end
function mu=sutherland(T)
mu=1.458e-6*T.^1.5/(T+110.4);
end
