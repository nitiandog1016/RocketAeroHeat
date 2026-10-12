function testEulerBoundaryReconstruction()
%TESTEULERBOUNDARYRECONSTRUCTION 独立仿射镜像场及总焓正值限制检查。
% 仿射夹具满足平面滑移镜像关系；不是实际曲面物理精度认证。
[x,r]=ndgrid(linspace(0,1,9),linspace(.1,1,7));
g=struct('X',x+.2*r,'R',r+.1*x,'nx',8,'nr',6);
rec=prepareEulerReconstruction(g,struct('boundaryReconstruction','slipGhost'));geo=rec.geometry;
[faces,source]=reconstructEulerPrimitive(affine(geo.centerX,geo.centerR),rec,false);
xx={geo.faceAX(1:end-1,:),geo.faceAX(2:end,:),geo.faceBX(:,1:end-1),geo.faceBX(:,2:end)};
rr={geo.faceAR(1:end-1,:),geo.faceAR(2:end,:),geo.faceBR(:,1:end-1),geo.faceBR(:,2:end)};
for k=1:4,assert(max(abs(faces{k}-affine(xx{k},rr{k})),[],'all')<1e-11);end
expected=affine(geo.areaCenterX,geo.areaCenterR);assert(max(abs(source-expected(:,:,4)),[],'all')<1e-11);
normal=[-.1,1]/sqrt(1.01);
assert(max(abs(faces{3}(:,1,2)*normal(1)+faces{3}(:,1,3)*normal(2)))<1e-11);
assert(nnz(rec.ghost)==g.nx && strcmp(rec.boundaryMethod,'slipGhost'));

%% 1. 轴线虚单元与常量状态：轴线零面积仍须有正确镜像位置
theta=linspace(pi/2,0,17)';contour=[-1,.1;.1*(cos(theta)-1),.1*sin(theta)];
g=buildNormalOffsetEulerGrid(contour,struct('axialCells',24,'radialCells',16));
rec=prepareEulerReconstruction(g,struct('boundaryReconstruction','slipGhost'));
assert(nnz(rec.ghost)==g.nx+g.nr);
q=repmat(reshape([1 0 0 1],1,1,4),g.nx,g.nr,1);
[faces,source]=reconstructEulerPrimitive(q,rec,true);
for k=1:4,assert(isequal(faces{k},q));end
assert(all(source==1,'all'));

%% 2. 接近动能主导的场：逐面保持常总焓和正压力，不改变单元原状态
[i,j]=ndgrid(1:g.nx,1:g.nr);phase=.7*i+1.3*j;rho=1+.1*sin(i);
u=4*cos(phase);v=4*sin(phase);H=8.00001;gamma=1.4;
p=rho.*(H-.5*(u.^2+v.^2))*(gamma-1)/gamma;q=cat(3,rho,u,v,p);before=q;
[faces,source,d]=reconstructEulerEnthalpy(q,rec,gamma);
assert(isequal(q,before) && d.enthalpyPositivityLimitedCellCount>0 && all(source>0,'all'));
for k=1:4
 f=faces{k};h=gamma/(gamma-1)*f(:,:,4)./f(:,:,1)+.5*(f(:,:,2).^2+f(:,:,3).^2);
 assert(max(abs(h-H),[],'all')<1e-12 && all(f(:,:,[1 4])>0,'all'));
end
bad=q;bad(:,:,4)=-1;reject(@()reconstructEulerEnthalpy(bad,rec,gamma),'RocketAeroHeat:EulerReconstruction:State');
bad=q;bad(:,:,4)=realmin;reject(@()reconstructEulerEnthalpy(bad,rec,gamma),'RocketAeroHeat:EulerReconstruction:Resolution');
fprintf('Boundary reconstruction passed: affine slip mirror, axis symmetry, positive static enthalpy and constant H.\n');
end
function q=affine(x,r)
% 切向坐标沿壁面，法向坐标在壁面为零；标量只随切向变化。
t=[1,.1]/sqrt(1.01);n=[-.1,1]/sqrt(1.01);s=t(1)*x+t(2)*r;
d=n(1)*x+n(2)*r-.098/sqrt(1.01);ut=20+.2*s;un=.3*d;
q=cat(3,2+.1*s,ut*t(1)+un*n(1),ut*t(2)+un*n(2),10+.2*s);
end
function reject(action,id)
try,action();error('Test:NotRejected','未拒绝非法数据');
catch problem,assert(strcmp(problem.identifier,id),problem.message);end
end
