function testEulerLinearReconstruction()
%TESTEULERLINEARRECONSTRUCTION 用独立锥流、均匀流及时间步检验重构求解器。
% 不把圆锥基准当作钝鼻或真实箭体热流精度证明。
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);ref=solveConicalShock(fs,10);
g=buildAxisymmetricEulerGrid([-1 0;0 0],struct('axialCells',24,'radialCells',16));
uniform=solveAxisymmetricEuler(g,fs,struct('spatialOrder',2));
assert(uniform.available && uniform.residual<1e-9 && uniform.spatialOrder==2);
assert(max(abs(uniform.pressure/fs.pressure-1),[],'all')<1e-12);
assert(uniform.conservation.identityError<1e-12);
errors=zeros(2,1);traceErrors=errors;
for k=1:2
 g=buildAxisymmetricEulerGrid([-1 tand(10);0 0], ...
     struct('axialCells',48*k,'radialCells',36*k,'radialStretch',1));
 s=solveAxisymmetricEuler(g,fs,struct('spatialOrder',2,'fluxScheme','HLL','CFL',.12,'maxIterations',20000));
 selected=g.wallX<-.35 & g.wallX>-.8;e=extractAxisymmetricEulerEdge(s);
 difference=s.wall.pressure(selected)/ref.pressure-1;
 % 使用逐点 RMS；有符号均值可能因正负误差抵消而随网格非单调。
 errors(k)=sqrt(mean(difference.^2));
 traceErrors(k)=sqrt(mean((s.wall.trace.pressure(selected)/ref.pressure-1).^2));
 assert(s.available && s.conservation.identityError<1e-10);
 assert(s.conservation.wallMassLeak<1e-12 && s.conservation.wallEnergyLeak<1e-12);
 assert(all(s.density>0,'all') && all(s.pressure>0,'all'));
 assert(~s.edgeStateAvailable && ~s.diagnostics.gridAccuracyVerified && ~e.available);
 assert(e.sample.fromReconstructedTrace && errors(k)<.015);
 fprintf('Linear cone %dx%d: pressure RMS=%g mean=%g trace RMS=%g residual=%g iterations=%d\n',g.nx,g.nr,errors(k),mean(difference),traceErrors(k),s.residual,s.iterations);
 if k==1,coarse=s;end
end
assert(errors(2)<errors(1) && traceErrors(2)<traceErrors(1));
% 由同网格解以另一伪时间步继续，仍要求完整重算残差并达到原阈值。
slow=solveAxisymmetricEuler(coarse.grid,fs,struct('spatialOrder',2,'fluxScheme','HLL', ...
    'CFL',.06,'initialState',coarse,'residualTolerance',1e-5,'maxIterations',20000));
assert(slow.available && slow.restart.sourceSpatialOrder==2);
assert(max(abs(slow.pressure-coarse.pressure),[],'all')/fs.pressure<1e-4);
mixed=solveAxisymmetricEuler(coarse.grid,fs,struct('spatialOrder',2,'fluxScheme','HLLenthalpy', ...
    'initialState',coarse,'maxIterations',20000));
assert(mixed.available && strcmp(mixed.restart.sourceFluxScheme,'HLL'));
assert(~mixed.edgeStateAvailable && mixed.conservation.identityError<1e-10);
fprintf('Limited linear Euler passed: uniform balance, cone refinement, conservative walls and pseudo-time restart.\n');
end
