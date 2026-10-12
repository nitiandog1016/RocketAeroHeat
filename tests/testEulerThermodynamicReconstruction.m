function testEulerThermodynamicReconstruction()
%TESTEULERTHERMODYNAMICRECONSTRUCTION 独立锥流检验镜像/总焓重构组合。
% 总焓偏差、温度、压力、速度和网格误差分别检查，不能互相替代。
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);ref=solveConicalShock(fs,10);
options=struct('spatialOrder',2,'boundaryReconstruction','slipGhost', ...
 'reconstructionVariables','enthalpy','fluxScheme','HLLenthalpy','maxIterations',16000);
g=buildAxisymmetricEulerGrid([-1 0;0 0],struct('axialCells',24,'radialCells',16));
uniform=solveAxisymmetricEuler(g,fs,options);
assert(uniform.available && uniform.residual<1e-9 && uniform.conservation.identityError<1e-12);
errors=zeros(2,3);saved=cell(2,1);
for k=1:2
 g=buildAxisymmetricEulerGrid([-1 tand(10);0 0],struct('axialCells',48*k,'radialCells',36*k,'radialStretch',1));
 s=solveAxisymmetricEuler(g,fs,options);e=extractAxisymmetricEulerEdge(s);selected=e.x<-.35 & e.x>-.8;
 errors(k,:)=[sqrt(mean((e.pressure(selected)/ref.pressure-1).^2)), ...
  sqrt(mean((e.temperature(selected)/ref.temperature-1).^2)),sqrt(mean((e.signedSpeed(selected)/ref.speed-1).^2))];
 assert(s.available && s.diagnostics.maximumTotalEnthalpyDeviation<1e-6);
 assert(s.conservation.identityError<1e-10 && s.conservation.wallMassLeak<1e-12 && s.conservation.wallEnergyLeak<1e-12);
 assert(all(errors(k,:)<[.015 .01 .002]));
 assert(~s.edgeStateAvailable && ~e.available && ~e.gridAccuracyVerified);
 assert(strcmp(s.boundaryReconstruction,'slipGhost') && strcmp(s.reconstructionVariables,'enthalpy'));
 assert(contains(s.method,'镜像') && contains(s.method,'rho/u/v/H'));
 fprintf('Thermodynamic cone %d: p RMS=%g T RMS=%g U RMS=%g H=%g\n',48*k,errors(k,:),s.diagnostics.maximumTotalEnthalpyDeviation);
 saved{k}=s;
end
assert(all(errors(2,:)<errors(1,:)));
comparison=compareAxisymmetricEulerEdges(saved{1},saved{2});assert(comparison.available && ~comparison.edgeStateAvailable);
bad=saved{2};bad.boundaryReconstruction='oneSided';reject(@()compareAxisymmetricEulerEdges(saved{1},bad));
bad=saved{2};bad.reconstructionVariables='primitive';reject(@()compareAxisymmetricEulerEdges(saved{1},bad));
% 切换重构变量可以作为同网格初场，但必须记录来源并重算收敛。
changed=options;changed.reconstructionVariables='primitive';changed.initialState=saved{1};changed.maxIterations=1;
restart=solveAxisymmetricEuler(saved{1}.grid,fs,changed);
assert(strcmp(restart.restart.sourceReconstructionVariables,'enthalpy') && ~restart.edgeStateAvailable);
fprintf('Thermodynamic Euler passed: independent cone state, grid refinement, H transport and method provenance.\n');
end
function reject(action)
try,action();error('Test:NotRejected','未拒绝算法不同的对比');
catch problem,assert(strcmp(problem.identifier,'RocketAeroHeat:EulerCompare:Reconstruction'),problem.message);end
end
