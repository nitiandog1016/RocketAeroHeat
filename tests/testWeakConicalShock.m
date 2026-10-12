function testWeakConicalShock()
%TESTWEAKCONICALSHOCK 独立定比热极限、实际弱锥流及半球/柱交界回归。
% 模块一：读取 SciPy 两变量无量纲解，不用被测函数产生期望值。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
ref=jsondecode(fileread(fullfile(root,'docs','validation','weak_cone_reference.json')));
fs=struct('Mach',ref.Mach,'temperature',ref.temperature,'R',ref.R, ...
    'gamma',ref.gamma,'pressure',1000);
fs.density=fs.pressure/(fs.R*fs.temperature);
fs.speed=fs.Mach*sqrt(fs.gamma*fs.R*fs.temperature);
actual=zeros(size(ref.rows));
for k=1:size(ref.rows,1)
    s=solveThermallyPerfectConicalShock(fs,ref.rows(k,1));
    assert(s.attached && s.available,'冻结比热极限的 %.3g 度锥未求解',ref.rows(k,1));
    actual(k,:)=[s.coneHalfAngleDeg,s.shockAngleDeg,s.pressureRatio,s.temperature];
end
error=abs(actual-ref.rows);
assert(max(error(:,2))<2e-6,'激波角与独立解不符');
assert(max(error(:,3))<2e-7,'压力比与独立解不符');
assert(max(error(:,4))<2e-5,'温度与独立解不符');
fprintf('Independent DOP853 weak cone comparison: max errors beta %.3g deg, p ratio %.3g, T %.3g K.\n',max(error(:,2:4)));

%% 模块二：真实 M4/20km 工况的小角度弱支、守恒和单调性
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=0;
fs=calcFreestream(cfg.flight);lastPressure=fs.pressure;
for angle=[.5 1 2 3 3.75 4 4.5 5 10]
    s=solveThermallyPerfectConicalShock(fs,angle);
    assert(s.attached && s.available,'实际工况 %.3g 度锥未求解',angle);
    assert(s.pressure>lastPressure);lastPressure=s.pressure;
    assert(abs(s.polarVelocityAtWall)<1e-3);
    assert(abs(s.energyResidual)<1e-7 && abs(s.shockMassResidual)<1e-9);
    assert(abs(s.shockMomentumResidual)<1e-9 && s.entropyIncrease>0);
end

%% 模块三：整模型走内部链路，交界带不得再因小角度求根失败而断开
[mesh,meta]=buildCanonicalBody('hemisphereCylinder',struct('resolution',48));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
cylinder=strcmp(meta.faceRegion,'cylinder');
% 第 32 环带的全部 96 面在旧版无外缘；用位置选取避免硬编码拓扑编号。
join=mesh.center(:,1)>0.4 & mesh.center(:,1)<0.4131;
assert(nnz(join)==96 && all(mesh.boundaryLayerEdge.algebraicallyValid(join)));
assert(all(mesh.surfaceStreamlines.connectedFace(join)));
assert(nnz(mesh.surfaceStreamlines.connectedFace & cylinder)>0.95*nnz(cylinder));
covered=sum(mesh.area(cylinder & mesh.heating.eligibleFace))/sum(mesh.area(cylinder));
assert(covered>0.85,'柱段热流仍有大范围缺失：%.3f',covered);
base=strcmp(meta.faceRegion,'base');
assert(all(isnan(mesh.heating.heatFlux(base))),'底面缺模型时不得填零');
fprintf('Hemisphere-cylinder join: all 96 faces have connected edge flow; cylinder computed heat area %.2f%%.\n',100*covered);
end
