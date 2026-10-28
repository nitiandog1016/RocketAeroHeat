function testAxisymmetricEntropyShockTable()
%TESTAXISYMMETRICENTROPYSHOCKTABLE 近似激波角/质量表的独立解析和收敛校核。
% 模块作用：检查 CR-112180 鼻端启动、几何/质量量纲、直锥解析不变量，
% 同时检查辅助压力与物理压力分开、热物性一致、域外拒绝和数值收敛。
% 本测试不证明真实激波形状、物理边界层闭合或实验热流精度。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
fs=makeState(220,5000,4);geometry=struct('noseRadius',.02,'baseRadius',.2, ...
    'coneAngleDeg',20,'tipX',1);
cfg.heating.airThermalModel='caloricallyPerfect';
cfg.externalFlow.entropyShock.stationCount=321;
timer=tic;table=solveAxisymmetricEntropyShockTable(fs,geometry,cfg);
assert(table.available,table.status);
assert(numel(table.massCoordinate)==321 && all(diff(table.surfaceArcLength)>0));
assert(all(diff(table.massCoordinate)>0) && all(diff(table.shockCrossingAngleDeg)<=0));
assert(table.massCoordinate(1)==0 && table.shockCrossingAngleDeg(1)==90);
assert(~table.geometryVerified && table.analyticGeometrySpecified && ...
    ~table.shockSurfaceAvailable && ~table.defaultSolverEnabled);
assert(all(strcmp(struct2cellKinds(table),'computed')));

%% 1. 独立完全气体正激波与修正牛顿压力，不复用激波核作为参考
g=fs.gamma;M=fs.Mach;Q=fs.density*fs.speed^2;
r=(g+1)*M^2/((g-1)*M^2+2);
p2=fs.pressure*(1+2*g/(g+1)*(M^2-1));T2=fs.temperature*(p2/fs.pressure)/r;
T0=fs.temperature*(1+(g-1)*M^2/2);
ps=p2*(T0/T2)^(g/(g-1));D=2*Q/(g+1);
A=((p2-fs.pressure)/Q)/(1+2*D/Q);
assert(abs(table.normalShockStaticPressure/p2-1)<1e-12);
assert(abs(table.stagnationPressure/ps-1)<1e-12 && ps>p2);
assert(abs(table.noseExpansion.A/A-1)<1e-12);
assert(abs(table.noseExpansion.cosAngleGradient*geometry.noseRadius/sqrt(2*A)-1)<1e-12);
t=table.surfaceArcLength/geometry.noseRadius;join=pi/2-deg2rad(geometry.coneAngleDeg);
mu=cos(min(t,join));expectedR=geometry.noseRadius*sin(min(t,join));
cone=t>join;mu(cone)=sind(geometry.coneAngleDeg);
expectedR(cone)=geometry.noseRadius*cosd(geometry.coneAngleDeg)+ ...
    (table.surfaceArcLength(cone)-geometry.noseRadius*join)*sind(geometry.coneAngleDeg);
assert(max(abs(table.surfaceRadius-expectedR))<1e-14);
assert(max(abs(table.physicalPressure-(fs.pressure+(ps-fs.pressure)*mu.^2)))<1e-9);
assert(max(abs(table.auxiliaryShockODEPressure-(fs.pressure+(p2-fs.pressure)*mu.^2)))<1e-9);
assert(table.physicalPressure(1)==ps && table.auxiliaryShockODEPressure(1)==p2);
assert(max(abs(table.physicalPressure-table.auxiliaryShockODEPressure))>1000);

%% 2. 鼻端展开、质量积分导数和直锥独立解析不变量
u=sind(table.shockCrossingAngleDeg);psi=table.massCoordinate/(fs.density*fs.speed*geometry.noseRadius^2);
q=t(2)^2;
assert(abs((1-u(2))/q/A-1)<1e-6);
assert(abs(psi(2)-(q/2-(A+1/6)*q^2/4))<1e-15);
expectedDerivative=fs.density*fs.speed*table.surfaceRadius.*u;
assert(max(abs(table.massCoordinateDerivative-expectedDerivative))<1e-13);
% 独立在物理S坐标对rho_inf*U_inf*r(S)*sinGamma作自适应积分，
% 不以ODE内部的Psi状态或差分自身作为质量守恒参考。
massIntegral=integral(@(s)massIntegrand(s,table,fs),0,table.surfaceArcLength(end), ...
    'Waypoints',table.surfaceArcLength(2:end-1),'RelTol',1e-8,'AbsTol',1e-12);
massIntegralError=abs(massIntegral/table.massCoordinate(end)-1);
assert(massIntegralError<1e-5);
% 辅助压力为常数的直锥上，(sin^2Gamma-u_inf^2)*Psi^(2D/Q) 守恒。
% 这是由完全气体 Psh 的闭式法向跳跃独立积分得到的精确ODE不变量。
uInf2=sind(geometry.coneAngleDeg)^2+cosd(geometry.coneAngleDeg)^2/M^2;
invariant=(u(cone).^2-uInf2).*psi(cone).^(2*D/Q);
assert(max(abs(invariant/invariant(1)-1))<5e-5);
far=geometry;far.baseRadius=1;farTable=solveAxisymmetricEntropyShockTable(fs,far,cfg);
assert(farTable.available,farTable.status);
% 有限长锥仍有可计算的熵层尾差，不能把有限端点当作无穷远。
farGap=sind(farTable.shockCrossingAngleDeg(end))^2-uInf2;
assert(farGap>0 && farGap<.01*(u(end)^2-uInf2));
farPsi=farTable.massCoordinate(end)/(fs.density*fs.speed*geometry.noseRadius^2);
assert(abs(farGap*farPsi^(2*D/Q)/invariant(1)-1)<5e-5);
fprintf('Caloric table: %d stations, terminal Gamma %.10f deg, cone invariant %.3g, %.2f s.\n', ...
    numel(t),table.shockCrossingAngleDeg(end),max(abs(invariant/invariant(1)-1)),toc(timer));
fprintf('Independent physical mass integral: relative error %.3g.\n',massIntegralError);

%% 3. 几何尺度协变及启动长度/压力表/积分容差收敛
scaled=geometry;scaled.noseRadius=3*geometry.noseRadius;scaled.baseRadius=3*geometry.baseRadius;
scaled.tipX=7;
scaledTable=solveAxisymmetricEntropyShockTable(fs,scaled,cfg);assert(scaledTable.available);
assert(max(abs(scaledTable.shockCrossingAngleDeg-table.shockCrossingAngleDeg))<1e-10);
assert(max(abs(scaledTable.massCoordinate/9-table.massCoordinate))<1e-12);
assert(max(abs(scaledTable.streamlineMetric/3-table.streamlineMetric))<1e-14);
fine=cfg;fine.externalFlow.entropyShock.noseStartRatio=5e-4;
fine.externalFlow.entropyShock.pressureTableNodes=513;
fine.externalFlow.entropyShock.relativeTolerance=1e-10;
fine.externalFlow.entropyShock.absoluteTolerance=1e-13;
fine.externalFlow.entropyShock.stationCount=321;
finer=solveAxisymmetricEntropyShockTable(fs,geometry,fine);assert(finer.available,finer.status);
compareTable(table,finer,2e-5);

%% 4. NASA7 同物性前向跃迁与导数/表格收敛，不把推广称原报告复现
thermal=cfg;thermal.heating.airThermalModel='thermallyPerfect';
thermal.externalFlow.entropyShock.pressureTableNodes=257;
thermalTable=solveAxisymmetricEntropyShockTable(fs,geometry,thermal);
assert(thermalTable.available,thermalTable.status);
assert(strcmp(thermalTable.thermalPropertyMode,'thermallyPerfect'));
assert(thermalTable.noseExpansion.pressureDerivativeRelativeError<1e-4);
jump=recoverShockEntropyEdgeState(fs,thermalTable.shockCrossingAngleDeg,fs.pressure,thermal);
assert(all(jump.shock.valid));
assert(max(abs(thermalTable.shockStaticPressure./jump.shock.pressure-1))<1e-5);
thermalFine=fine;thermalFine.heating.airThermalModel='thermallyPerfect';
thermalFine.externalFlow.entropyShock.pressureTableNodes=1025;
thermalFiner=solveAxisymmetricEntropyShockTable(fs,geometry,thermalFine);
assert(thermalFiner.available,thermalFiner.status);compareTable(thermalTable,thermalFiner,5e-5);
fprintf('NASA7 table: terminal Gamma %.10f deg, pressure derivative uncertainty %.3g.\n', ...
    thermalTable.shockCrossingAngleDeg(end),thermalTable.noseExpansion.pressureDerivativeRelativeError);

%% 5. 真实不适用条件保留拒绝，不回退到其他压力或补质量
oblique=fs;oblique.alpha=1;refused=solveAxisymmetricEntropyShockTable(oblique,geometry,cfg);
assert(~refused.available && isempty(refused.massCoordinate));
low=makeState(220,5000,2);refused=solveAxisymmetricEntropyShockTable(low,geometry,cfg);
assert(~refused.available && isempty(refused.massCoordinate));
outside=geometry;outside.coneAngleDeg=10;failed=false;
try,solveAxisymmetricEntropyShockTable(fs,outside,cfg);catch problem,failed=contains(problem.identifier,'Geometry');end
assert(failed);
save(fullfile(root,'outputs','v52_validation','shock_table_results.mat'), ...
    'table','finer','thermalTable','thermalFiner','fs','geometry');
fprintf('V52 AXISYMMETRIC SHOCK TABLE PASSED: source startup, caloric invariant, scale and numerical convergence; no physical shock-surface or experimental claim.\n');
end

function fs=makeState(T,p,M)
R=287.05;g=1.4;fs=struct('temperature',T,'pressure',p,'R',R,'gamma',g, ...
    'Mach',M,'density',p/(R*T),'speed',M*sqrt(g*R*T),'alpha',0,'beta',0);
end
function kinds=struct2cellKinds(table)
names=fieldnames(table.fieldMetadata);kinds=cell(numel(names),1);
for j=1:numel(names),kinds{j}=table.fieldMetadata.(names{j}).valueKind;end
end
function compareTable(coarse,fine,tolerance)
% 比较相同物理站点；细表有密集输出，仅比较远离启动的共同区间。
% 两次启动改变了鼻部对数输出网格。直锥部分使用相同的物理站，
% 不将跨稀疏网格插值误差误算为ODE的启动/容差收敛误差。
query=coarse.surfaceArcLength;
join=coarse.geometry.noseRadius*(pi/2-deg2rad(coarse.geometry.coneAngleDeg));
keep=query>=join;
angle=interp1(fine.surfaceArcLength,fine.shockCrossingAngleDeg,query(keep),'pchip');
mass=interp1(fine.surfaceArcLength,fine.massCoordinate,query(keep),'pchip');
fprintf('ODE/start/pressure-grid convergence: angle %.3g deg, relative mass %.3g.\n', ...
    max(abs(angle-coarse.shockCrossingAngleDeg(keep))),max(abs(mass./coarse.massCoordinate(keep)-1)));
assert(max(abs(angle-coarse.shockCrossingAngleDeg(keep)))<tolerance);
assert(max(abs(mass./coarse.massCoordinate(keep)-1))<tolerance);
end
function value=massIntegrand(s,table,fs)
Rn=table.geometry.noseRadius;alpha=table.geometry.coneAngleDeg;
join=Rn*(pi/2-deg2rad(alpha));r=Rn*sin(min(s/Rn,join/Rn));
cone=s>join;r(cone)=Rn*cosd(alpha)+(s(cone)-join)*sind(alpha);
u=interp1(table.surfaceArcLength,sind(table.shockCrossingAngleDeg),s,'pchip');
value=fs.density*fs.speed*r.*u;
end
