function reconstruction = prepareEulerReconstruction(grid,options)
%PREPAREEULERRECONSTRUCTION 预计算真实坐标的加权最小二乘梯度系数。
% 默认只使用域内邻居；可选滑移镜像虚单元，用物理边界关系辅助梯度。
% 系数每张网格只算一次。顺序为 i-、i+、j-、j+。
if nargin<2,options=struct;end
boundary='oneSided';if isfield(options,'boundaryReconstruction'),boundary=options.boundaryReconstruction;end
boundary=validatestring(boundary,{'oneSided','slipGhost'});
geo=axisymmetricCellGeometry(grid.X,grid.R);nx=grid.nx;nr=grid.nr;
ids=reshape(1:nx*nr,nx,nr);neighbors=cat(3,ids([1 1:end-1],:), ...
    ids([2:end end],:),ids(:,[1 1:end-1]),ids(:,[2:end end]));
neighbors=reshape(neighbors,[],4);center=ids(:);valid=neighbors~=center;
x=geo.centerX(:);r=geo.centerR(:);
dx=x(neighbors)-x;dr=r(neighbors)-r;
ghost=false(size(dx));normalX=zeros(size(dx));normalR=zeros(size(dx));
if strcmp(boundary,'slipGhost')
    % 镜像单元位于物理壁面/轴线另一侧；只用于梯度，不新增控制体。
    % 面法向由节点几何计算，因此零旋转面积的轴线也有有限法向。
    % 标量/切向速度偶延拓，法向速度奇延拓；这是平面滑移边界近似，
    % 不代表已经解析真实曲壁的法向压力梯度或壁面涡量。
    wx=diff(grid.X(:,1));wr=diff(grid.R(:,1));wl=hypot(wx,wr);
    cells=ids(:,1);normalX(cells,3)=-wr./wl;normalR(cells,3)=wx./wl;
    distance=(geo.faceBX(:,1)-x(cells)).*normalX(cells,3)+(geo.faceBR(:,1)-r(cells)).*normalR(cells,3);
    dx(cells,3)=2*distance.*normalX(cells,3);dr(cells,3)=2*distance.*normalR(cells,3);ghost(cells,3)=true;
    if isfield(grid,'boundaryTopology') && strcmp(grid.boundaryTopology,'axisAtNose')
        cells=ids(end,:);normalR(cells,2)=1;dx(cells,2)=0;dr(cells,2)=-2*r(cells);ghost(cells,2)=true;
    end
    valid=valid | ghost;
end
w=zeros(size(dx));w(valid)=1./(dx(valid).^2+dr(valid).^2);
a=sum(w.*dx.^2,2);b=sum(w.*dx.*dr,2);d=sum(w.*dr.^2,2);determinant=a.*d-b.^2;
if any(determinant<=1e-12)
    error('RocketAeroHeat:EulerReconstruction:Stencil','重构邻域退化，无法拟合二维梯度。');
end
reconstruction.geometry=geo;reconstruction.neighbors=neighbors;reconstruction.interiorFace=neighbors~=center;
reconstruction.stencilActive=valid;
reconstruction.ghost=ghost;reconstruction.ghostNormalX=normalX;reconstruction.ghostNormalR=normalR;
reconstruction.boundaryMethod=boundary;
reconstruction.coefficientX=w.*(d.*dx-b.*dr)./determinant;
reconstruction.coefficientR=w.*(a.*dr-b.*dx)./determinant;
fx=cat(3,geo.faceAX(1:end-1,:),geo.faceAX(2:end,:),geo.faceBX(:,1:end-1),geo.faceBX(:,2:end));
fr=cat(3,geo.faceAR(1:end-1,:),geo.faceAR(2:end,:),geo.faceBR(:,1:end-1),geo.faceBR(:,2:end));
reconstruction.deltaX=reshape(fx,[],4)-x;reconstruction.deltaR=reshape(fr,[],4)-r;
reconstruction.sourceDeltaX=geo.areaCenterX(:)-x;reconstruction.sourceDeltaR=geo.areaCenterR(:)-r;
reconstruction.nx=nx;reconstruction.nr=nr;
reconstruction.smoothingLengthCubed=(sqrt(geo.area(:))/(max(grid.X,[],'all')-min(grid.X,[],'all'))).^3;
reconstruction.method='物理坐标反距离平方加权最小二乘；域内四邻居，边界单边拟合';
if strcmp(boundary,'slipGhost')
    reconstruction.method='物理坐标反距离平方加权最小二乘；域内邻居与壁面/轴线平面滑移镜像虚单元';
end
end
