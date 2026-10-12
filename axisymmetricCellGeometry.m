function geometry = axisymmetricCellGeometry(X,R)
%AXISYMMETRICCELLGEOMETRY 一般四边形的轴对称积分矩与面取样位置。
% 模块用途：线性重构必须以体积加权中心为基准，径向压力源使用面积
% 加权中心；两者通常不同。所有积分省略共同的 2*pi。
% 每个四边形分成两个三角形，用二次精确三点公式计算 r、x*r、r^2；
% X 在每个单元内平移，避免大坐标相减造成多边形矩的消减误差。
validateattributes(X,{'numeric'},{'real','finite','2d'});
validateattributes(R,{'numeric'},{'real','finite','nonnegative','size',size(X)});
if any(size(X)<2),error('RocketAeroHeat:EulerGeometry:Size','至少需要 2×2 个节点。');end
xa=X(1:end-1,1:end-1);ra=R(1:end-1,1:end-1);
A=cat(3,zeros(size(xa)),ra);
B=cat(3,X(2:end,1:end-1)-xa,R(2:end,1:end-1));
C=cat(3,X(2:end,2:end)-xa,R(2:end,2:end));
D=cat(3,X(1:end-1,2:end)-xa,R(1:end-1,2:end));
[a1,v1,x1,r1,ax1,ar1]=moments(A,B,C);
[a2,v2,x2,r2,ax2,ar2]=moments(A,C,D);
area=a1+a2;volume=v1+v2;
if any(a1<=0,'all') || any(a2<=0,'all') || any(volume<=0,'all')
    error('RocketAeroHeat:EulerGeometry:Orientation','单元三角形面积或旋转体积非正。');
end
geometry.volume=volume;geometry.area=area;
geometry.centerX=xa+(x1+x2)./volume;geometry.centerR=(r1+r2)./volume;
geometry.areaCenterX=xa+(ax1+ax2)./area;geometry.areaCenterR=(ar1+ar2)./area;
% 直线面的 r 加权中心：参数 t 的精确均值为 (r1+2*r2)/(3*(r1+r2))。
% 轴线面的旋转面积为零，取几何中点作为有限的重构位置。
[geometry.faceAX,geometry.faceAR]=faceCenter(X(:,1:end-1),R(:,1:end-1),X(:,2:end),R(:,2:end));
[geometry.faceBX,geometry.faceBR,geometry.faceBFraction]=faceCenter(X(1:end-1,:),R(1:end-1,:),X(2:end,:),R(2:end,:));
geometry.method='三角形二次精确积分；旋转体积中心、面积中心和 r 加权面中心分别计算';
end
function [area,volume,mx,mr,ax,ar]=moments(A,B,C)
u=B-A;v=C-A;area=.5*(u(:,:,1).*v(:,:,2)-u(:,:,2).*v(:,:,1));
volume=zeros(size(area));mx=volume;mr=volume;
points=cat(4,(4*A+B+C)/6,(A+4*B+C)/6,(A+B+4*C)/6);
for k=1:3
    x=points(:,:,1,k);r=points(:,:,2,k);w=area/3;
    volume=volume+w.*r;mx=mx+w.*x.*r;mr=mr+w.*r.^2;
end
ax=area.*(A(:,:,1)+B(:,:,1)+C(:,:,1))/3;
ar=area.*(A(:,:,2)+B(:,:,2)+C(:,:,2))/3;
end
function [x,r,t]=faceCenter(x1,r1,x2,r2)
t=.5*ones(size(r1));active=(r1+r2)>0;
t(active)=(r1(active)+2*r2(active))./(3*(r1(active)+r2(active)));
x=x1+t.*(x2-x1);r=r1+t.*(r2-r1);
end
