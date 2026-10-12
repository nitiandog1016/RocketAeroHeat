function grid = buildNormalOffsetEulerGrid(contour,options)
%BUILDNORMALOFFSETEULERGRID 沿母线节点法向展开的钝鼻 Euler 试验网格。
% 模块用途：减轻 X=常数网格在圆钝鼻附近的严重倾斜。
% 输入轮廓按 X 递增，鼻端为 R=0；末端切线须近似垂直主轴。
% i 从尾部走向鼻端，j 从壁面走向远场；i 的最后一条边是鼻前轴线，
% 不是入口。远场包括鼻前圆弧状边界，要求求解器支持一般面法向。
% 保留全部轮廓折点，不光滑 STL；凹角可能导致法线相交，必须拒绝坏网格。
if nargin<2,options=struct;end
validateattributes(contour,{'numeric'},{'real','finite','2d','ncols',2});
if size(contour,1)<3 || any(diff(contour(:,1))<=0) || any(contour(1:end-1,2)<=0) || contour(end,2)~=0
    error('RocketAeroHeat:NormalGrid:Contour','要求正半径连续母线及 +X 端唯一的零半径鼻尖。');
end
L=contour(end,1)-contour(1,1);
nx=option(options,'axialCells',120);nr=option(options,'radialCells',64);
outer=option(options,'outerRadius',max(.65*L,3*max(contour(:,2))));
stretch=option(options,'radialStretch',1.8);
blend=option(options,'normalBlendLength',inf);
refinement=option(options,'wallRefinement',1);
validateattributes(nx,{'numeric'},{'scalar','integer','>=',8,'<=',2000});
validateattributes(nr,{'numeric'},{'scalar','integer','>=',8,'<=',1000});
validateattributes(outer,{'numeric'},{'real','finite','scalar','>',max(contour(:,2))});
validateattributes(stretch,{'numeric'},{'real','finite','scalar','>=',1,'<=',3});
validateattributes(blend,{'numeric'},{'real','scalar','positive'});
validateattributes(refinement,{'numeric'},{'scalar','integer','>=',1,'<=',8});
last=contour(end,:)-contour(end-1,:);
if last(1)/norm(last)>.1
    error('RocketAeroHeat:NormalGrid:Nose','法向网格当前要求圆钝鼻端；不适用于尖锥。');
end

%% 1. 按弧长补充壁面节点，在原始折线上建立节点平均法向
arc=[0;cumsum(hypot(diff(contour(:,1)),diff(contour(:,2))))];
subdivided=arc(1:end-1)+diff(arc)*(0:refinement-1)/refinement;
locations=unique([arc;subdivided(:);linspace(0,arc(end),nx+1)']);
locations=locations([true;diff(locations)>1e-12*L]);
wall=interp1(arc,contour,locations,'linear');
tangent=[wall(2,:)-wall(1,:);wall(3:end,:)-wall(1:end-2,:);wall(end,:)-wall(end-1,:)];
tangent=tangent./vecnorm(tangent,2,2);
normal=[-tangent(:,2),tangent(:,1)];
if isfinite(blend)
    % 可选鼻部过渡：仅在前端逐渐采用法向，其余区域保持竖直网格线，
    % 避免将肩部凹角法线一直延长到远场。只改变网格线，不平滑壁面。
    distance=arc(end)-locations;weight=max(0,1-(distance/blend).^2).^2;
    nxNormal=cummax(max(0,normal(:,1).*weight));
    normal=[nxNormal,sqrt(max(0,1-nxNormal.^2))];
end
% 鼻端节点位于对称轴，外侧节点也必须留在轴上；不改变任何单元速度。
normal(end,:)=[1,0];
eta=linspace(0,1,nr+1).^stretch;
offset=(outer-max(contour(:,2)))*eta;
X=wall(:,1)+normal(:,1)*offset;R=wall(:,2)+normal(:,2)*offset;
R(end,:)=0;

%% 2. 对每个一般四边形精确积分面积及旋转体积（省略共同 2*pi）
xa=X(1:end-1,1:end-1);ra=R(1:end-1,1:end-1);
xb=X(2:end,1:end-1);rb=R(2:end,1:end-1);
xc=X(2:end,2:end);rc=R(2:end,2:end);
xd=X(1:end-1,2:end);rd=R(1:end-1,2:end);
ab=xa.*rb-xb.*ra;bc=xb.*rc-xc.*rb;cd=xc.*rd-xd.*rc;da=xd.*ra-xa.*rd;
area=.5*(ab+bc+cd+da);
volume=((ra+rb).*ab+(rb+rc).*bc+(rc+rd).*cd+(rd+ra).*da)/6;
% 两个网格方向的有向面面积，均指向对应索引增大的方向。
meanA=.5*(R(:,1:end-1)+R(:,2:end));
Ax=diff(R,1,2).*meanA;Ar=-diff(X,1,2).*meanA;
meanB=.5*(R(1:end-1,:)+R(2:end,:));
Bx=-diff(R,1,1).*meanB;Br=diff(X,1,1).*meanB;
metric=max([abs(diff(Ax,1,1)+diff(Bx,1,2)), ...
    abs(diff(Ar,1,1)+diff(Br,1,2)-area)],[],'all')/max(area,[],'all');
% 正面积不足以排除凹四边形；各相邻边叉积也必须为正。
c1=(xb-xa).*(rd-ra)-(rb-ra).*(xd-xa);
c2=(xc-xb).*(ra-rb)-(rc-rb).*(xa-xb);
c3=(xd-xc).*(rb-rc)-(rd-rc).*(xb-xc);
c4=(xa-xd).*(rc-rd)-(ra-rd).*(xc-xd);
if any(R<0,'all') || any(volume<=0,'all') || any(min(cat(3,c1,c2,c3,c4),[],3)<=0,'all') || metric>1e-9
    error('RocketAeroHeat:NormalGrid:Intersection', ...
        '法向偏移产生相交/凹单元或度量误差；请使用原 X 截面网格，不能直接用于该外形。');
end
grid=struct('X',X,'R',R,'volume',volume,'meridionalArea',area, ...
    'axialArea',Ax,'axialAreaR',Ar,'radialAreaX',Bx,'radialAreaR',Br, ...
    'centerX',.25*(xa+xb+xc+xd),'centerR',.25*(ra+rb+rc+rd), ...
    'wallX',.5*(wall(1:end-1,1)+wall(2:end,1)), ...
    'wallR',.5*(wall(1:end-1,2)+wall(2:end,2)), ...
    'bodyWall',true(size(wall,1)-1,1),'nx',size(wall,1)-1,'nr',nr, ...
    'contour',contour,'length',L,'outerRadius',max(R,[],'all'), ...
    'radialStretch',stretch,'metricError',metric,'boundaryTopology','axisAtNose');
grid.method='母线节点法向偏移；一般四边形精确轴对称度量；鼻前轴线对称边界（试验）';
grid.normalBlendLength=blend;
grid.wallRefinement=refinement;
if isfinite(blend)
    grid.method=sprintf('鼻部法向—竖直过渡网格，弧长范围 %.4g m；保留壁面折点；鼻前轴线（试验）',blend);
end
grid.method=sprintf('%s；母线每段细分 %d',grid.method,refinement);
end
function value=option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
