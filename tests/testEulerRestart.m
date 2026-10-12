function testEulerRestart()
%TESTEULERRESTART 验证重启、通量切换、稳态一致和过期/非法状态拒绝。
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);
g=buildAxisymmetricEulerGrid([-1 tand(10);0 0],struct('axialCells',32,'radialCells',24));
base=solveAxisymmetricEuler(g,fs,struct('fluxScheme','HLL'));
same=solveAxisymmetricEuler(g,fs,struct('fluxScheme','HLL','initialState',base));
assert(same.converged && same.iterations==1 && same.restart.used);
assert(max(abs(same.pressure-base.pressure),[],'all')/fs.pressure<1e-12);
assert(strcmp(same.restart.sourceFluxScheme,'HLL'));
changed=solveAxisymmetricEuler(g,fs,struct('initialState',base,'maxIterations',1));
assert(~changed.converged && changed.restart.sourceConverged);
continued=solveAxisymmetricEuler(g,fs,struct('initialState',changed));
cold=solveAxisymmetricEuler(g,fs);
assert(continued.available && cold.available && ~continued.restart.sourceConverged);
assert(max(abs(continued.pressure-cold.pressure),[],'all')/fs.pressure<1e-4);
assert(continued.diagnostics.maximumTotalEnthalpyDeviation<1e-6);
assert(continued.diagnostics.maximumTotalEnthalpyDeviation<base.diagnostics.maximumTotalEnthalpyDeviation/1000);
assert(contains(continued.method,'上风总焓') && ~continued.edgeStateAvailable);
newGrid=buildAxisymmetricEulerGrid([-1 tand(10);0 0],struct('axialCells',34,'radialCells',24));
mustReject(@()solveAxisymmetricEuler(newGrid,fs,struct('initialState',base)),'RestartGrid');
flight=cfg.flight;flight.Mach=flight.Mach+1;newFs=calcFreestream(flight);
mustReject(@()solveAxisymmetricEuler(g,newFs,struct('initialState',base)),'RestartFlight');
bad=base;bad.freestream.alpha=5;
mustReject(@()solveAxisymmetricEuler(g,fs,struct('initialState',bad)),'RestartFlight');
bad=base;bad.density(1)=-1;
mustReject(@()solveAxisymmetricEuler(g,fs,struct('initialState',bad)),'RestartState');
fprintf('Euler restart passed: cold/warm equivalence, changed flux re-evaluation, enthalpy transport and stale state rejection.\n');
end
function mustReject(f,id)
try,f();error('Test:NotRejected','非法输入未拒绝');
catch problem,assert(strcmp(problem.identifier,['RocketAeroHeat:Euler:' id]));end
end
