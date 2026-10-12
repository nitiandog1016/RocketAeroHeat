function testNormalOffsetEulerGrid()
%TESTNORMALOFFSETEULERGRID 一般四边形的独立体积积分、几何守恒及轴线边界。
theta=linspace(pi/2,0,33)';c=[-1,.1;.1*(cos(theta)-1),.1*sin(theta)];
g=buildNormalOffsetEulerGrid(c,struct('axialCells',48,'radialCells',16));
assert(g.metricError<1e-10 && all(g.volume>0,'all') && all(g.R(end,:)==0));
assert(all(g.axialArea(end,:)==0) && all(g.axialAreaR(end,:)==0));
assert(all(g.bodyWall));
% 每个四边形分成两个三角形，以三角形面积乘重心半径独立积分体积。
A=cat(3,g.X(1:end-1,1:end-1),g.R(1:end-1,1:end-1));
B=cat(3,g.X(2:end,1:end-1),g.R(2:end,1:end-1));
C=cat(3,g.X(2:end,2:end),g.R(2:end,2:end));
D=cat(3,g.X(1:end-1,2:end),g.R(1:end-1,2:end));
v=triangle(A,B,C)+triangle(A,C,D);
assert(max(abs(v-g.volume)./v,[],'all')<1e-10);
% 静止常压 Euler 场：一般法向压力通量与径向几何源精确抵消。
p=3.7;gam=1.4;Q=reshape([1,0,0,p/(gam-1)],1,1,4);
fx=axisymmetricEulerFlux(repmat(Q,g.nx+1,g.nr),repmat(Q,g.nx+1,g.nr), ...
    g.axialArea,g.axialAreaR,gam,'HLL');
fb=axisymmetricEulerFlux(repmat(Q,g.nx,g.nr+1),repmat(Q,g.nx,g.nr+1), ...
    g.radialAreaX,g.radialAreaR,gam,'HLL');
r=diff(fx,1,1)+diff(fb,1,2);r(:,:,3)=r(:,:,3)-p*g.meridionalArea;
assert(max(abs(r),[],'all')<1e-12);
% 每条原始轮廓折点均保留，不能通过平滑几何获得好看的结果。
for j=1:size(c,1),assert(min(hypot(g.X(:,1)-c(j,1),g.R(:,1)-c(j,2)))<1e-12);end
blended=buildNormalOffsetEulerGrid(c,struct('normalBlendLength',.1));
assert(all(blended.volume>0,'all') && blended.metricError<1e-10);
assert(contains(blended.method,'过渡') && blended.normalBlendLength==.1);
assert(all(blended.X(blended.wallX<-.3,1)==blended.X(blended.wallX<-.3,end)));
fine=buildNormalOffsetEulerGrid(c,struct('axialCells',48,'radialCells',16,'wallRefinement',2));
assert(isequal(fine.contour,g.contour) && fine.nx>g.nx && fine.wallRefinement==2);
assert(abs(hypot(fine.wallX(end),fine.wallR(end))/hypot(g.wallX(end),g.wallR(end))-.5)<1e-10);
old=buildAxisymmetricEulerGrid(c);new=buildAxisymmetricEulerGrid(c,struct('wallRefinement',2));
assert(isequal(old.contour,new.contour) && new.nx>old.nx);
assert(abs(new.wallX(find(new.bodyWall,1,'last'))/old.wallX(find(old.bodyWall,1,'last'))-.5)<1e-10);
reject(@()buildNormalOffsetEulerGrid([-1 .2;-.5 .1;0 0]),'RocketAeroHeat:NormalGrid:Nose');
reject(@()buildNormalOffsetEulerGrid([-1 .1;-.5 0;0 0]),'RocketAeroHeat:NormalGrid:Contour');
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);
bad=g;bad.R(end,1)=.01;
reject(@()solveAxisymmetricEuler(bad,fs),'RocketAeroHeat:Euler:AxisBoundary');
fprintf('Normal-offset Euler grid passed: independent volumes, constant-pressure balance, axis and contour checks.\n');
end
function volume=triangle(A,B,C)
u=B-A;v=C-A;area=.5*(u(:,:,1).*v(:,:,2)-u(:,:,2).*v(:,:,1));
volume=area.*(A(:,:,2)+B(:,:,2)+C(:,:,2))/3;
end
function reject(f,id)
try,f();error('Test:NotRejected','未拒绝非法网格');
catch problem,assert(strcmp(problem.identifier,id),problem.message);end
end
