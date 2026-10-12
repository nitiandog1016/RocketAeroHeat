function testBowShockPreview()
% Billig 球形鼻端预览：顶点、脱体距离、曲率、网格与适用域。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();');
fs=calcFreestream(cfg.flight);
R=0.05; tip=[1 0.2 -0.1];
m.V=[tip;0.8 0.2 -0.1;0.9 0.24 -0.1];
m.stagnation.noseRadiusEstimate=struct('available',true,'radius',R);
m.stagnationShock.available=true;
b=estimateBowShockSurface(m,fs);
expected=R*0.143*exp(3.24/fs.Mach^2);
assert(b.available && abs(b.standOffDistance-expected)<1e-12);
assert(abs(b.apexPosition(1)-tip(1)-expected)<1e-12);
assert(all(abs(b.apexPosition(2:3)-tip(2:3))<1e-12));
assert(size(b.V,2)==3 && size(b.F,2)==3);
assert(all(b.F(:)>=1 & b.F(:)<=size(b.V,1)));
assert(all(diff(b.meridianX)<0));
bodyX=tip(1)-R+sqrt(max(0,R^2-b.meridianRadius.^2));
assert(all(b.meridianX>bodyX));
assert(contains(b.source,'10.2514/3.28969'));
slow=fs; slow.Mach=2.5;
slowPreview=estimateBowShockSurface(m,slow);
assert(~slowPreview.available);
tilt=fs; tilt.alpha=10;
tiltedPreview=estimateBowShockSurface(m,tilt);
assert(~tiltedPreview.available);
fprintf('Bow-shock preview passed: Delta/Rn=%.4f, Rc/Rn=%.4f.\n', ...
    b.standOffDistance/R,b.shockCurvatureRadius/R);
end
