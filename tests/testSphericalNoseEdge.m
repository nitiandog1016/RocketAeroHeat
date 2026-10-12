function testSphericalNoseEdge()
%TESTSPHERICALNOSEEDGE 球鼻外缘的独立解析、守恒、几何域及物性拒绝校核。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=20;
cfg.flight.beta=0;cfg.flight.Twall=300;
fs=calcFreestream(cfg.flight);
[mesh,meta]=buildCanonicalBody('hemisphereCylinder',struct('resolution',24));
mesh=prepare(mesh,cfg); before=mesh;
edge=solveSphericalNoseEdge(mesh,fs,cfg);
assert(isequaln(mesh,before),'独立模型不得改写网格或旧求解结果');
valid=edge.algebraicallyValid;
assert(edge.available && nnz(valid)>30);
assert(all(valid([848 850 852])),'必须覆盖已经定位的球鼻外流缺失带');
assert(~any(valid(strcmp(meta.faceRegion,'cylinder') | strcmp(meta.faceRegion,'base'))));
assert(all(edge.polarAngleDeg(valid)<=60+1e-10));
assert(all(isnan(edge.pressure(~valid))) && all(isnan(edge.temperature(~valid))));
assert(max(abs(edge.pressure(valid)-edge.density(valid).*fs.R.*edge.temperature(valid)))<1e-8);
assert(max(abs(edge.energyResidual(valid)))<1e-10);
assert(max(abs(edge.entropyResidual(valid)))<1e-8);
assert(max(abs(sum(edge.velocityVector(valid,:).*mesh.normal(valid,:),2)))<1e-9);
assert(max(abs(vecnorm(edge.velocityVector(valid,:),2,2)-edge.speed(valid)))<1e-9);
[~,order]=sort(edge.polarAngleDeg(valid));p=edge.pressure(valid);U=edge.speed(valid);
assert(all(diff(p(order))<=1e-8) && all(diff(U(order))>=-1e-8));
% 将来流准确对准一个球鼻面元的径向方向，检查真实零速度极限。
r=mesh.center-edge.sphereCenter;r=r./vecnorm(r,2,2);
id=find(edge.geometryFace & r(:,1)>cosd(20),1);
aligned=fs;aligned.unitFlowDirection=-r(id,:);
stagnant=solveSphericalNoseEdge(mesh,aligned,cfg);
assert(stagnant.algebraicallyValid(id) && stagnant.speed(id)==0);
assert(stagnant.pressure(id)==stagnant.totalPressure(id));
assert(stagnant.temperature(id)==stagnant.totalTemperature(id));
assert(all(stagnant.velocityVector(id,:)==0));

%% 2. 定比热独立公式校核；不调用被测模块内部的正激波函数作参考
caloric=cfg;caloric.heating.airThermalModel='caloricallyPerfect';
c=solveSphericalNoseEdge(mesh,fs,caloric);v=c.algebraicallyValid;
gamma=fs.gamma;M=fs.Mach;R=fs.R;
p2=fs.pressure*(2*gamma*M^2-(gamma-1))/(gamma+1);
M2sq=((gamma-1)*M^2+2)/(2*gamma*M^2-(gamma-1));
p02=p2*(1+(gamma-1)*M2sq/2)^(gamma/(gamma-1));
T0=fs.temperature+(fs.speed^2)/(2*gamma*R/(gamma-1));
center=[max(mesh.V(:,1))-.1,0,0];r=mesh.center(v,:)-center;
mu=(r*(-fs.unitFlowDirection(:)))./vecnorm(r,2,2);
expectedP=fs.pressure+(p02-fs.pressure)*mu.^2;
expectedT=T0*(expectedP/p02).^((gamma-1)/gamma);
expectedU=sqrt(2*gamma*R/(gamma-1)*(T0-expectedT));
assert(max(abs(c.pressure(v)./expectedP-1))<1e-12);
assert(max(abs(c.temperature(v)./expectedT-1))<1e-12);
assert(max(abs(c.speed(v)./expectedU-1))<1e-12);

%% 3. 几何变换不改场量；非球鼻、特征面、拓扑问题及过大攻角拒绝
shift=mesh;shift.V=2*mesh.V+[7 -3 5];shift=prepare(shift,cfg);
translated=solveSphericalNoseEdge(shift,fs,cfg);
assert(isequal(valid,translated.algebraicallyValid));
assert(max(abs(translated.pressure(valid)./edge.pressure(valid)-1))<1e-9);
bad=mesh;bad.V(:,2)=1.2*bad.V(:,2);bad=prepare(bad,cfg);
bad.stagnation=struct('noseRadiusEstimate',edge.geometryEstimate);
rejected=solveSphericalNoseEdge(bad,fs,cfg);assert(~rejected.available);
marked=mesh;marked.localGeometry.featureFace(848)=true;
marked.localGeometry.topologyIssue(850)=true;
protected=solveSphericalNoseEdge(marked,fs,cfg);
assert(~protected.algebraicallyValid(848) && ~protected.algebraicallyValid(850));
flight=cfg.flight;flight.alpha=35;
rejected=solveSphericalNoseEdge(mesh,calcFreestream(flight),cfg);assert(~rejected.available);
[coarse,~]=buildCanonicalBody('bluntedCone',struct('resolution',24));
coarse=prepare(coarse,cfg);rejected=solveSphericalNoseEdge(coarse,fs,cfg);
assert(~rejected.available,'粗钝锥不能因为给定模型半径就通过球面识别');

%% 4. 热物性越界不得隐式回退；范围收紧必须有明确输出
flight=cfg.flight;flight.Mach=10;
hot=solveSphericalNoseEdge(mesh,calcFreestream(flight),cfg);
assert(~hot.available && all(isnan(hot.pressure)));
narrow=cfg;narrow.externalFlow.sphericalNoseMaxPolarAngleDeg=30;
limited=solveSphericalNoseEdge(mesh,fs,narrow);
assert(nnz(limited.algebraicallyValid)<nnz(valid));
assert(~any(limited.algebraicallyValid & ~valid));
fprintf('Spherical nose edge passed: %d faces; missing-band faces restored locally; analytic perfect-gas, energy/entropy, geometry and scope checks passed.\n',nnz(valid));
fprintf('This is a standalone numerical/analytic check, not interface or experimental validation.\n');
end
function mesh=prepare(mesh,cfg)
% 从真实顶点重新生成几何；不伪造法向或球面识别结果。
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
end
