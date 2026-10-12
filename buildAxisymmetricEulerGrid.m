function grid = buildAxisymmetricEulerGrid(contour,options)
%BUILDAXISYMMETRICEULERGRID 构建含轴对称面积权重的贴体四边形网格。
% contour=[X,R]，按 X 递增，末端为 +X 鼻尖且 R=0；单位 m。
% 只处理单值回转体外形，尾端直接截断，不包含底部和尾迹。
% 所有体积、面积均省略共同的 2*pi，守恒方程必须使用相同约定。
if nargin<2, options=struct; end
validateattributes(contour,{'numeric'},{'real','finite','2d','ncols',2});
if size(contour,1)<2 || any(diff(contour(:,1))<=0) || ...
        any(contour(:,2)<0) || contour(end,2)~=0
    error('RocketAeroHeat:Euler:Contour','轮廓 X 必须严格递增，半径非负，鼻尖半径必须为零。');
end
L=contour(end,1)-contour(1,1);
nx=option(options,'axialCells',120); nr=option(options,'radialCells',64);
outer=option(options,'outerRadius',max(0.65*L,3*max(contour(:,2))));
upstream=option(options,'upstreamLength',0.15*L);
stretch=option(options,'radialStretch',1.8);
refinement=option(options,'wallRefinement',1);
validateattributes(nx,{'numeric'},{'scalar','integer','>=',8,'<=',2000});
validateattributes(nr,{'numeric'},{'scalar','integer','>=',8,'<=',1000});
validateattributes(outer,{'numeric'},{'real','finite','scalar','>',max(contour(:,2))});
validateattributes(upstream,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(stretch,{'numeric'},{'real','finite','scalar','>=',1,'<=',3});
validateattributes(refinement,{'numeric'},{'scalar','integer','>=',1,'<=',8});

%% 1. 轴向和径向节点：鼻尖处强制有网格线，保持轴/壁边界分界
% 分别划分机体与鼻前区，避免鼻尖和浮点近重复节点产生极薄网格。
% 保留所有母线折点，不把肩部等几何特征在粗网格上平滑掉。
nb=max(2,round(nx*L/(L+upstream)));nf=max(2,nx-nb);
% 鼻前轴向单元也需解析脱体激波间距；仅径向加密会留下粗大鼻前单元。
% 首格不大于鼻端最小母线间距，向远场缓慢增长，增长比不超过 1.15。
firstStep=min([L/nx,upstream/nf,min(diff(contour(:,1)))/refinement]);
nf=max(nf,ceil(log(1+.15*upstream/firstStep)/log(1.15)));
weights=1.15.^(0:nf-1);noseAhead=upstream*[0,cumsum(weights)/sum(weights)];
subdivided=contour(1:end-1,1)+diff(contour(:,1))*(0:refinement-1)/refinement;
x=unique([linspace(contour(1,1),contour(end,1),nb+1), ...
    contour(end,1)+noseAhead,contour(:,1).',subdivided(:).']).';
% 只合并舍入误差尺度的邻近点，优先保留用户给定的几何折点。
for k=numel(x)-1:-1:1
    if x(k+1)-x(k)<1e-12*L
        if any(contour(:,1)==x(k)),x(k+1)=[];else,x(k)=[];end
    end
end
rb=interp1(contour(:,1),contour(:,2),min(x,contour(end,1)),'linear');
eta=linspace(0,1,nr+1).^stretch;
R=rb+(outer-rb).*eta; X=repmat(x,1,nr+1);
dx=diff(x); loL=R(1:end-1,1:end-1); loR=R(2:end,1:end-1);
hiL=R(1:end-1,2:end); hiR=R(2:end,2:end);

%% 2. 精确积分线性边界的旋转体积与加权面法向
grid.volume=dx/6.*(hiL.^2+hiL.*hiR+hiR.^2-loL.^2-loL.*loR-loR.^2);
grid.meridionalArea=dx/2.*(hiL+hiR-loL-loR);
grid.axialArea=0.5*diff(R.^2,1,2);
rMean=0.5*(R(1:end-1,:)+R(2:end,:));
grid.radialAreaX=-diff(R,1,1).*rMean;
grid.radialAreaR=dx.*rMean;
grid.X=X; grid.R=R;
grid.centerX=0.5*(x(1:end-1)+x(2:end))*ones(1,nr);
grid.centerR=0.25*(loL+loR+hiL+hiR);
grid.wallX=0.5*(x(1:end-1)+x(2:end));
grid.wallR=0.5*(rb(1:end-1)+rb(2:end));
grid.bodyWall=grid.wallX<contour(end,1);
grid.contour=contour; grid.length=L; grid.outerRadius=outer;
grid.nx=numel(x)-1; grid.nr=nr;
grid.radialStretch=stretch;
grid.wallRefinement=refinement;
% 恒压的轴向净通量为零，径向净通量必须等于几何源项。
grid.metricError=max([abs(diff(grid.axialArea,1,1)+diff(grid.radialAreaX,1,2)), ...
    abs(diff(grid.radialAreaR,1,2)-grid.meridionalArea)],[],'all')/max(grid.meridionalArea,[],'all');
if any(grid.volume<=0,'all') || grid.metricError>1e-10
    error('RocketAeroHeat:Euler:Grid','网格体积或几何守恒检查失败。');
end
grid.method=sprintf('分片线性贴体四边形，精确轴对称面积/体积积分（省略共同 2π）；母线每段细分 %d',refinement);
end
function value=option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
