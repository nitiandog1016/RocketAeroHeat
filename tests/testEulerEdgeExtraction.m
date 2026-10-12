function testEulerEdgeExtraction()
%TESTEULEREDGEEXTRACTION 独立解析场验证弧长、投影、非均匀导数及网格对比。
% 人造解析场仅检验提取算子；不作为 Euler PDE 收敛或物理精度的证据。
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);
c=[-1 tand(10);-.83 .83*tand(10);-.43 .43*tand(10);0 0];
a=fixture(c,24,fs);e=extractAxisymmetricEulerEdge(a);
assert(all(diff(e.s)>0) && all(diff(e.x)<0) && e.s(1)>0);
assert(max(abs(e.signedSpeed-(100+20*e.s+5*e.s.^2)))<1e-10);
assert(max(abs(e.velocityGradient-(20+10*e.s)))<1e-8);
assert(max(abs(e.pressureGradient-fs.pressure*(.1+.1*e.s)))<1e-7);
assert(max(abs(e.pressure./(e.density*fs.R)-e.temperature))<1e-10);
assert(max(abs(e.sample.totalEnthalpy-e.totalEnthalpy-.5*e.sample.normalVelocity.^2))<1e-8);
assert(max(abs(e.velocityX.*e.normal(:,1)+e.velocityR.*e.normal(:,2)))<1e-10);
assert(~e.available && ~e.gridAccuracyVerified && ~any(e.flags.feature));
assert(contains(e.method,'第一层') && contains(e.gradientMethod,'非均匀'));
% 一阶取样的法向分量、壁压差必须可见，不能因投影而消失。
assert(all(e.diagnostics.normalVelocityFraction>.2));
assert(max(abs(e.diagnostics.wallPressureRelativeDifference-.03))<1e-12);
unfinished=a;unfinished.available=false;unfinished.converged=false;
unfinished.wall.pressure(e.wallIndex(1))=-1;
u=extractAxisymmetricEulerEdge(unfinished);
assert(~u.inputConverged && ~u.available && u.flags.wallFluxPressureNonpositive(1));
assert(u.wallFluxPressure(1)==-1 && u.pressure(1)>0);

%% 重构壁面迹：使用 r 加权位置，独立多项式验证取样和非均匀导数
traceCase=a;traceCase.spatialOrder=2;g=a.grid;
dx=diff(g.X(:,1));dr=diff(g.R(:,1));lengths=hypot(dx,dr);
r0=g.R(1:end-1,1);r1=g.R(2:end,1);fraction=.5*ones(g.nx,1);
active=r0+r1>0;fraction(active)=(r0(active)+2*r1(active))./(3*(r0(active)+r1(active)));
t=struct('x',g.X(1:end-1,1)+fraction.*dx,'radius',r0+fraction.*dr, ...
    'fractionFromTail',fraction,'density',a.density(:,1),'pressure',a.pressure(:,1), ...
    'velocityX',a.velocityX(:,1),'velocityR',a.velocityR(:,1));
ids=e.wallIndex;ss=cumsum(lengths(ids))-fraction(ids).*lengths(ids);
speed=100+20*ss+5*ss.^2;vn=.25*speed;
t.velocityX(ids)=speed.*e.tangent(:,1)+vn.*e.normal(:,1);
t.velocityR(ids)=speed.*e.tangent(:,2)+vn.*e.normal(:,2);
t.pressure(ids)=fs.pressure*(1+.1*ss+.05*ss.^2);
traceCase.wall.trace=t;et=extractAxisymmetricEulerEdge(traceCase);
assert(et.sample.fromReconstructedTrace && contains(et.method,'重构'));
assert(abs(et.s(1)/lengths(ids(1))-2/3)<1e-12);
assert(max(abs(et.signedSpeed-speed))<1e-10);
assert(max(abs(et.velocityGradient-(20+10*ss)))<1e-8);
assert(max(abs(et.pressureGradient-fs.pressure*(.1+.1*ss)))<1e-7);
assert(max(abs(et.sample.normalVelocity-vn))<1e-10 && ~et.available);
bad=traceCase;bad.wall.trace.pressure(1)=0;
reject(@()extractAxisymmetricEulerEdge(bad),'RocketAeroHeat:EulerEdge:Trace');

%% 1. 真实几何折角及逆流处分段，不能跨越不适用区域差分
kink=fixture([-1 .1;-.5 .1;0 0],48,fs);k=extractAxisymmetricEulerEdge(kink);
assert(nnz(k.flags.feature)==2 && all(isnan(k.velocityGradient(k.flags.feature))));
assert(max(abs(k.velocityGradient(k.flags.gradientAvailable)-(20+10*k.s(k.flags.gradientAvailable))))<1e-8);
reverse=a;i=e.wallIndex(round(numel(e.s)/2));
reverse.velocityX(i,1)=-a.velocityX(i,1);reverse.velocityR(i,1)=-a.velocityR(i,1);
r=extractAxisymmetricEulerEdge(reverse);
assert(nnz(r.flags.reverseFlow)==1 && all(isnan(r.velocityGradient(r.flags.reverseFlow))));
bad=a;bad.density(1)=0;reject(@()extractAxisymmetricEulerEdge(bad),'RocketAeroHeat:EulerEdge:State');
bad=a;bad.freestream.alpha=3;reject(@()extractAxisymmetricEulerEdge(bad),'RocketAeroHeat:EulerEdge:Incidence');

%% 2. 同外形双网格对比；解析线性场能精确插值，不把重复计算当加密
a=linearFixture(a);b=linearFixture(fixture(c,48,fs));
comparison=compareAxisymmetricEulerEdges(a,b);
assert(comparison.available && ~comparison.gridAccuracyVerified && ~comparison.edgeStateAvailable);
assert(comparison.axialNodesChanged);
assert(comparison.signedSpeed.maximumAbsoluteDifference<1e-10);
assert(comparison.pressure.maximumAbsoluteDifference<1e-9);
assert(comparison.velocityGradient.maximumAbsoluteDifference<1e-8);
reject(@()compareAxisymmetricEulerEdges(a,a),'RocketAeroHeat:EulerCompare:SameGrid');
bad=b;bad.available=false;reject(@()compareAxisymmetricEulerEdges(a,bad),'RocketAeroHeat:EulerCompare:Unconverged');
bad=b;bad.freestream.Mach=6;reject(@()compareAxisymmetricEulerEdges(a,bad),'RocketAeroHeat:EulerCompare:Flight');
bad=b;bad.fluxScheme='HLL';reject(@()compareAxisymmetricEulerEdges(a,bad),'RocketAeroHeat:EulerCompare:Scheme');
bad=b;bad.spatialOrder=2;reject(@()compareAxisymmetricEulerEdges(a,bad),'RocketAeroHeat:EulerCompare:Reconstruction');
bad=fixture([-1 .3;0 0],48,fs);reject(@()compareAxisymmetricEulerEdges(a,bad),'RocketAeroHeat:EulerCompare:Geometry');
fprintf('Euler edge extraction passed: analytic gradients, slip projection, feature gaps and grid comparison.\n');
end

function a=fixture(contour,n,fs)
% 独立指定解析近壁场；故意含法向速度，以检测是否隐藏非切向误差。
g=buildAxisymmetricEulerGrid(contour,struct('axialCells',n,'radialCells',16));
a=struct('grid',g,'freestream',fs,'density',ones(g.nx,g.nr)*fs.density, ...
    'pressure',ones(g.nx,g.nr)*fs.pressure,'velocityX',zeros(g.nx,g.nr), ...
    'velocityR',zeros(g.nx,g.nr),'thermodynamicModel','caloricallyPerfect', ...
    'method','解析测试夹具，非 PDE 解','fluxScheme','HLLenthalpy','available',true,'converged',true);
a.boundary.clear=true;a.wall.pressure=ones(g.nx,1)*fs.pressure;
ids=flipud(find(g.bodyWall & g.wallR>0));dx=diff(g.X(:,1));dr=diff(g.R(:,1));ds=hypot(dx,dr);
s=cumsum(ds(ids))-.5*ds(ids);tx=-dx(ids)./ds(ids);tr=-dr(ids)./ds(ids);
u=100+20*s+5*s.^2;v=.25*u;
a.velocityX(ids,1)=u.*tx+v.*tr;a.velocityR(ids,1)=u.*tr-v.*tx;
a.pressure(ids,1)=fs.pressure*(1+.1*s+.05*s.^2);a.wall.pressure(ids)=1.03*a.pressure(ids,1);
end
function a=linearFixture(a)
e=extractAxisymmetricEulerEdge(a);ids=e.wallIndex;s=e.s;
a.velocityX(ids,1)=(100+20*s).*e.tangent(:,1);
a.velocityR(ids,1)=(100+20*s).*e.tangent(:,2);
a.pressure(ids,1)=a.freestream.pressure*(1+.1*s);a.wall.pressure(ids)=a.pressure(ids,1);
end
function reject(action,id)
try,action();error('Test:NotRejected','未拒绝不兼容数据');
catch problem,assert(strcmp(problem.identifier,id),problem.message);end
end
