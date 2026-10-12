function testAerodynamicLoads()
% 解析方柱验证：压力抵消、力矩平移、参考面积缩放和正负攻角对称。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
m.V=[0 -1 -1;0 1 -1;0 1 1;0 -1 1;1 -1 -1;1 1 -1;1 1 1;1 -1 1];
m.F=[1 2 6;1 6 5;2 3 7;2 7 6;3 4 8;3 8 7;4 1 5;4 5 8;1 4 3;1 3 2;5 6 7;5 7 8];
evalc('m=calcFaceGeometry(m); m=fixNormalDirection(m,cfg); m=buildSurfaceCoordinates(m,cfg);');
m=analyzeLocalGeometry(m,cfg);
cfg.aerodynamics=struct('referenceArea',4,'referenceLength',1,'momentReference',[0 0 0]);
flight=cfg.flight; flight.alpha=0; flight.beta=0;
m=modifiedNewtonianPressure(m,calcFreestream(flight),cfg);
uniform=m; uniform.pressure.pressureIncrement(:)=100;
uniform=integrateAerodynamicLoads(uniform,cfg);
assert(norm(uniform.aerodynamicLoads.forceBody)<1e-10);
assert(norm(uniform.aerodynamicLoads.momentBody)<1e-10);
nose=m; nose.pressure.pressureIncrement(:)=0;
nose.pressure.pressureIncrement(nose.localGeometry.isNoseCap)=100;
nose=integrateAerodynamicLoads(nose,cfg); r=nose.aerodynamicLoads;
assert(norm(r.forceBody-[-400 0 0])<1e-10 && norm(r.momentBody)<1e-10);
shifted=cfg; shifted.aerodynamics.momentReference=[0 1 0];
ms=integrateAerodynamicLoads(nose,shifted);
assert(norm(ms.aerodynamicLoads.momentBody-[0 0 -400])<1e-10);
scaled=cfg; scaled.aerodynamics.referenceArea=8;
ms=integrateAerodynamicLoads(nose,scaled);
assert(norm(ms.aerodynamicLoads.forceBody-r.forceBody)==0);
assert(abs(ms.aerodynamicLoads.CD-r.CD/2)<1e-12);
forces=zeros(3,3); moments=zeros(3,3);
for k=1:3
    flight.alpha=(k-2)*20;
    mm=modifiedNewtonianPressure(m,calcFreestream(flight),cfg);
    mm=integrateAerodynamicLoads(mm,cfg);
    forces(k,:)=mm.aerodynamicLoads.forceBody;
    moments(k,:)=mm.aerodynamicLoads.momentBody;
end
assert(abs(forces(2,2))+abs(forces(2,3))<1e-10);
assert(abs(forces(1,1)-forces(3,1))<1e-10);
assert(abs(forces(1,3)+forces(3,3))<1e-10 && forces(3,3)>0);
assert(abs(moments(1,2)+moments(3,2))<1e-10);
bad=cfg; bad.aerodynamics.referenceArea=0;
mustFail(@()integrateAerodynamicLoads(m,bad));
bad=m; bad.F(1,:)=[];
mustFail(@()integrateAerodynamicLoads(bad,cfg));
bad=m; bad.pressure.geometry.center(1,1)=bad.pressure.geometry.center(1,1)+1;
mustFail(@()integrateAerodynamicLoads(bad,cfg));
fprintf('Aerodynamic loads tests passed: analytic force, moment shift, normalization, symmetry, invalid inputs.\n');
end
function mustFail(f)
try, f(); catch, return; end
error('预期错误未发生。');
end
