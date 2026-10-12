function testAxisymmetricEuler()
%TESTAXISYMMETRICEULER 验证均匀流、守恒、非法工况与独立锥流基准。
evalc('cfg=createCaseConfig();'); fs=calcFreestream(cfg.flight);
fs.alpha=0; fs.beta=0;
options=struct('axialCells',64,'radialCells',48,'radialStretch',1,'maxIterations',16000);
grid=buildAxisymmetricEulerGrid([-1 0;0 0],options);
uniform=solveAxisymmetricEuler(grid,fs,options);
assert(uniform.available && uniform.residual<1e-9);
assert(max(abs(uniform.pressure/fs.pressure-1),[],'all')<1e-12);
assert(max(abs(uniform.velocityR),[],'all')<1e-12);
assert(uniform.conservation.identityError<1e-12);
bad=fs;bad.alpha=1;
try, solveAxisymmetricEuler(grid,bad,options); error('Test:NotRejected','未拒绝攻角');
catch err,assert(strcmp(err.identifier,'RocketAeroHeat:Euler:Incidence'));end
grid=buildAxisymmetricEulerGrid([-1 tand(10);0 0],options);
unfinished=solveAxisymmetricEuler(grid,fs,struct('maxIterations',1));
assert(~unfinished.available && ~unfinished.converged && ~unfinished.edgeStateAvailable);
sol=solveAxisymmetricEuler(grid,fs,options);
ref=solveConicalShock(fs,10);
selected=grid.wallX<-.35 & grid.wallX>-.8;
err=mean(sol.wall.pressure(selected))/ref.pressure-1;
fprintf('Euler cone: available=%d, iterations=%d, residual=%.5g, wall error=%.5g, outer=%.5g\n', ...
    sol.available,sol.iterations,sol.residual,err,sol.boundary.outerPressureDisturbance);
assert(sol.available && abs(err)<0.03);
assert(sol.conservation.identityError<1e-11);
assert(sol.conservation.wallMassLeak<1e-12 && sol.conservation.wallEnergyLeak<1e-12);
assert(all(sol.density>0,'all') && all(sol.temperature>0,'all'));
assert(~sol.edgeStateAvailable);
% 单元压力（而非反射壁面通量压力）的网格误差应随加密下降。
coarse=solveAxisymmetricEuler(buildAxisymmetricEulerGrid([-1 tand(10);0 0], ...
    struct('axialCells',32,'radialCells',24,'radialStretch',1)),fs);
fine=solveAxisymmetricEuler(buildAxisymmetricEulerGrid([-1 tand(10);0 0], ...
    struct('axialCells',128,'radialCells',96,'radialStretch',1)),fs);
sequence={coarse,sol,fine};cellError=zeros(3,1);
for j=1:3
    s=sequence{j}; mask=s.grid.wallX<-.35 & s.grid.wallX>-.8;
    cellError(j)=abs(mean(s.wall.cellPressure(mask))/ref.pressure-1);
    assert(s.available);
end
assert(all(cellError(2:3)<0.7*cellError(1:2)));
[~,i]=min(abs(fine.grid.wallX+.7));[~,j]=max(-diff(fine.pressure(i,:)));
shockAngle=atand(mean(fine.grid.centerR(i,j:j+1))/(-fine.grid.wallX(i)));
assert(abs(shockAngle-ref.shockAngleDeg)<.8);
fprintf('Refinement cell pressure errors: %.5g %.5g %.5g; shock angle %.5g vs %.5g deg\n',cellError,shockAngle,ref.shockAngleDeg);
% 改变伪时间步不应改变收敛稳态；较近远场必须被诊断为不合格。
slow=solveAxisymmetricEuler(grid,fs,struct('CFL',.35));
assert(slow.available && max(abs(slow.pressure-sol.pressure),[],'all')/fs.pressure<1e-4);
near=solveAxisymmetricEuler(buildAxisymmetricEulerGrid([-1 tand(10);0 0], ...
    struct('axialCells',64,'radialCells',48,'outerRadius',.2)),fs);
assert(~near.boundary.clear && ~near.available);
fprintf('PASS axisymmetric Euler kernel and independent cone benchmark\n');
end
