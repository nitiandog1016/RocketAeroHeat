function mesh = integrateAerodynamicLoads(mesh,cfg)
%INTEGRATEAERODYNAMICLOADS 对闭合表面的压力载荷进行气动力与力矩积分。
% 输入：最终网格 mesh（含 localGeometry、pressure），可选 cfg。
% 输出：mesh.aerodynamicLoads；不修改网格、压力或局部几何。
% 单位：长度 m，面积 m²，力 N，力矩 N·m；系数无量纲。
%
% 面元力 = -(p-pinf)*A*n；面元力矩 = (面心-参考点) × 面元力。
% 法向须一致向外；仅压力贡献，不包含摩擦阻力、推力或重力。
% 体轴：+X 鼻向，+Y 横向，+Z 向上；力矩按右手定则。
% 风轴：阻力正向沿空气流动方向，零攻角升力 +Z、侧向力 +Y。
% 对一般攻角/侧滑，升力方向为 [sin(alpha),0,cos(alpha)]；
% 侧向方向由 cross(升力方向,飞行方向) 得到。
% 注意 [阻力,侧力,升力] 是投影顺序，不是右手坐标基。
%
% cfg.aerodynamics.referenceArea：参考面积，默认最大轴向截面面积。
% cfg.aerodynamics.referenceLength：参考长度，默认 X 方向全长。
% cfg.aerodynamics.momentReference：[x0 y0 z0]，默认尾端几何中心。
% 默认参考点不是质心！所有力矩系数统一用 q*Sref*Lref 归一化。

if nargin<2, cfg=struct; end
%% 1. 检查当前几何与压力是否相容，拒绝开口或绕序不一致的网格
if ~all(isfield(mesh,{'V','F','localGeometry','pressure'}))
    error('RocketAeroHeat:Loads:Input','请先完成局部几何和压力计算。');
end
V=double(mesh.V); F=double(mesh.F); g=mesh.localGeometry; p=mesh.pressure;
validateattributes(V,{'numeric'},{'real','finite','2d','ncols',3,'nonempty'});
validateattributes(F,{'numeric'},{'real','finite','integer','positive','2d','ncols',3,'nonempty','<=',size(V,1)});
nf=size(F,1);
validateattributes(p.pressureIncrement,{'numeric'},{'real','finite','size',[nf 1]});
validateattributes(p.dynamicPressure,{'numeric'},{'real','finite','scalar','positive'});
edges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
[~,~,group]=unique(sort(edges,2),'rows'); count=accumarray(group,1);
direction=2*double(edges(:,2)>edges(:,1))-1;
balance=accumarray(group,direction);
if any(count~=2) || any(balance~=0)
    error('RocketAeroHeat:Loads:Topology','整体载荷积分要求闭合、流形且绕序一致的网格。');
end
a=V(F(:,1),:); b=V(F(:,2),:)-a; c=V(F(:,3),:)-a;
raw=cross(b,c,2); area=vecnorm(raw,2,2)/2;
if any(area<=0), error('RocketAeroHeat:Loads:Geometry','存在退化面元。'); end
normal=raw./(2*area); center=a+(b+c)/3;
scale=max(max(V)-min(V)); tol=max(1e-10*scale,64*eps(max(abs(V(:)))));
if ~isequal(size(g.normal),[nf 3]) || ~isequal(size(g.center),[nf 3]) || ...
        ~isequal(size(g.area),[nf 1]) || any(vecnorm(g.normal-normal,2,2)>1e-10) || ...
        any(vecnorm(g.center-center,2,2)>tol) || any(abs(g.area-area)>1e-10*area)
    error('RocketAeroHeat:Loads:StaleGeometry','几何已改变，请重建局部几何和压力。');
end
if isfield(p,'geometry')
    if ~isequaln(p.geometry.normal,g.normal) || ~isequaln(p.geometry.center,g.center) || ...
            ~isequaln(p.geometry.area,g.area)
        error('RocketAeroHeat:Loads:StalePressure','压力对应旧网格，请重新计算压力。');
    end
end
volume=sum(sum((center-mean(V,1)).*normal,2).*area)/3;
if volume<=0, error('RocketAeroHeat:Loads:Orientation','外法向校验失败，请执行 fixNormalDirection。'); end

%% 2. 确定并记录参考量，不把几何中心自动当作质心
options=struct;
if isfield(cfg,'aerodynamics'), options=cfg.aerodynamics; end
if isfield(options,'referenceArea')
    S=options.referenceArea; areaSource='用户指定';
elseif isfield(mesh,'param') && isfield(mesh.param,'sectionAreaYZ')
    S=max(abs(mesh.param.sectionAreaYZ)); areaSource='最大轴向截面面积';
else
    % 没有截面数据时使用 YZ 投影凸包面积；该后备值明确标记来源。
    yz=unique(V(:,2:3),'rows'); hull=convhull(yz(:,1),yz(:,2));
    S=polyarea(yz(hull,1),yz(hull,2)); areaSource='YZ 投影凸包面积（后备值）';
end
L=max(V(:,1))-min(V(:,1)); lengthSource='模型轴向全长';
if isfield(options,'referenceLength'), L=options.referenceLength; lengthSource='用户指定'; end
aft=V(abs(V(:,1)-min(V(:,1)))<=tol,:);
origin=[min(V(:,1)) mean(aft(:,2:3),1)]; pointSource='尾端节点几何中心（非质心）';
% STL 端面可能略有倾斜，最小 X 节点的平均值未必在中心线上。
% 有截面重构信息时沿截面中心线外推到尾端，避免选中单个边缘节点。
if isfield(mesh,'param') && all(isfield(mesh.param,{'sectionX','sectionCenterYZ'})) ...
        && numel(mesh.param.sectionX)>=2
    yz=interp1(mesh.param.sectionX,mesh.param.sectionCenterYZ,min(V(:,1)),'linear','extrap');
    origin=[min(V(:,1)) yz]; pointSource='尾端截面中心线外推（非质心）';
end
if isfield(options,'momentReference'), origin=options.momentReference; pointSource='用户指定'; end
validateattributes(S,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(L,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(origin,{'numeric'},{'real','finite','vector','numel',3}); origin=origin(:)';

%% 3. 积分面元力与相对参考点的力矩
faceForce=-(p.pressureIncrement.*area).*normal;
faceMoment=cross(center-origin,faceForce,2);
force=sum(faceForce,1); moment=sum(faceMoment,1);

%% 4. 按项目攻角与侧滑定义建立风轴投影，并检查方向一致性
fs=p.freestream;
validateattributes(fs.alpha,{'numeric'},{'real','finite','scalar'});
validateattributes(fs.beta,{'numeric'},{'real','finite','scalar'});
flight=[cosd(fs.alpha)*cosd(fs.beta),sind(fs.beta),-sind(fs.alpha)*cosd(fs.beta)];
lift=[sind(fs.alpha),0,cosd(fs.alpha)];
side=cross(lift,flight); drag=-flight;
flow=fs.unitFlowDirection(:)';
if numel(flow)~=3 || any(~isfinite(flow)) || norm(flow)==0 || norm(drag-flow/norm(flow))>1e-10
    error('RocketAeroHeat:Loads:FlowConvention','攻角、侧滑角与空气方向不一致。');
end
wind=force*[drag;side;lift]'; denominator=p.dynamicPressure*S;

%% 5. 保存分量、系数、参考量和模型范围，供统一 GUI 读取
r=struct;
r.faceForce=faceForce; r.faceMoment=faceMoment;
r.forceBody=force; r.momentBody=moment;
r.forceWind=wind; r.drag=wind(1); r.sideForce=wind(2); r.lift=wind(3);
r.forceCoefficientBody=force/denominator;
r.momentCoefficientBody=moment/(denominator*L);
r.CD=wind(1)/denominator; r.CY=wind(2)/denominator; r.CL=wind(3)/denominator;
r.windDirections=struct('drag',drag,'side',side,'lift',lift);
r.reference=struct('area',S,'length',L,'momentPoint',origin, ...
    'areaSource',areaSource,'lengthSource',lengthSource,'pointSource',pointSource);
r.dynamicPressure=p.dynamicPressure;
r.description='闭合表面表压积分；仅压力贡献，不含摩擦阻力';
mesh.aerodynamicLoads=r;
end
