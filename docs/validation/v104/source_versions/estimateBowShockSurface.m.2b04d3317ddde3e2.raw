function preview = estimateBowShockSurface(mesh,freestream)
%ESTIMATEBOWSHOCKSURFACE Billig 球形钝鼻经验式的局部三维激波预览。
% 模块作用：对零攻角、球形钝鼻、3<=M<=17 的完全气体工况，
% 估算脱体距离与激波顶点曲率，并将 r<=Rn 的轴对称曲线旋转成曲面。
% Billig, F. S., "Shock-wave shapes around spherical- and cylindrical-nosed
% bodies", Journal of Spacecraft and Rockets 4(6), 822–823 (1967),
% DOI: 10.2514/3.28969。下式只选球形钝鼻分支：
% Delta/Rn=0.143*exp(3.24/M^2)
% Rc/Rn=1.143*exp(0.54/(M-1)^1.2)
% 远场斜率使用 Mach 角。曲面只截取鼻端 r<=Rn，不代表整枚火箭激波。

%% 1. 几何和工况适用性检查
preview=struct('available',false, ...
    'method','Billig 1967 球形钝鼻经验激波几何关联式', ...
    'scope','零攻角球形钝鼻、3<=M<=17、近似完全气体；仅 r<=Rn 的鼻端局部预览', ...
    'source','Billig (1967), J. Spacecraft Rockets 4(6), 822–823, DOI:10.2514/3.28969', ...
    'status','当前几何或工况不满足经验关联式', ...
    'standOffDistance',NaN,'shockCurvatureRadius',NaN, ...
    'noseRadius',NaN,'machAngleDeg',NaN, ...
    'V',zeros(0,3),'F',zeros(0,3));
if ~isfield(mesh,'stagnation') || ...
        ~isfield(mesh.stagnation,'noseRadiusEstimate') || ...
        ~mesh.stagnation.noseRadiusEstimate.available
    preview.status='STL 鼻端未通过球形钝鼻多尺度检查'; return;
end
if ~isfield(mesh,'stagnationShock') || ~mesh.stagnationShock.available
    preview.status='当前没有可用的驻点法向激波热力状态'; return;
end
if freestream.Mach<3 || freestream.Mach>17 || ...
        hypot(freestream.alpha,freestream.beta)>1e-8 || ...
        abs(freestream.gamma-1.4)>0.05
    preview.status='Billig 球形鼻端预览要求 M=3~17、零攻角及近似空气定比热';
    return;
end
R=mesh.stagnation.noseRadiusEstimate.radius;
if ~isfinite(R) || R<=0
    preview.status='球形鼻尖半径无效'; return;
end

%% 2. 球形钝鼻脱体距离、激波顶点曲率与轴对称母线
M=freestream.Mach;
delta=R*0.143*exp(3.24/M^2);
Rc=R*1.143*exp(0.54/(M-1)^1.2);
machAngle=asin(1/M);
tanAngle=tan(machAngle);
cotSquared=1/tanAngle^2;
Vbody=double(mesh.V);
xTip=max(Vbody(:,1));
tip=Vbody(abs(Vbody(:,1)-xTip)<=max(1e-9,1e-8*R),:);
axisYZ=mean(tip(:,2:3),1);
radialStations=25; circumferentialPoints=48;
r=linspace(0,R,radialStations+1)';
z=(r*tanAngle/Rc).^2;
% 用 z/(sqrt(1+z)+1) 避免鼻尖处 sqrt(1+z)-1 的消减误差。
x=xTip+delta-Rc*cotSquared*z./(sqrt(1+z)+1);

%% 3. 旋转母线并三角化；单个顶点封闭最前端，避免退化面
theta=(0:circumferentialPoints-1)*2*pi/circumferentialPoints;
V=zeros(1+radialStations*circumferentialPoints,3);
V(1,:)=[x(1),axisYZ];
for j=1:radialStations
    ids=2+(j-1)*circumferentialPoints:1+j*circumferentialPoints;
    V(ids,1)=x(j+1);
    V(ids,2)=axisYZ(1)+r(j+1)*cos(theta);
    V(ids,3)=axisYZ(2)+r(j+1)*sin(theta);
end
F=zeros(circumferentialPoints*(1+2*(radialStations-1)),3);
k=0;
for c=1:circumferentialPoints
    next=mod(c,circumferentialPoints)+1;
    k=k+1; F(k,:)=[1,1+c,1+next];
end
for j=1:radialStations-1
    a=2+(j-1)*circumferentialPoints;
    b=a+circumferentialPoints;
    for c=0:circumferentialPoints-1
        next=mod(c+1,circumferentialPoints);
        k=k+1; F(k,:)=[a+c,b+c,b+next];
        k=k+1; F(k,:)=[a+c,b+next,a+next];
    end
end
preview.available=true;
preview.status='球形钝鼻经验激波面已生成；仅供工程几何预览';
preview.standOffDistance=delta;
preview.shockCurvatureRadius=Rc;
preview.noseRadius=R;
preview.machAngleDeg=rad2deg(machAngle);
preview.tipPosition=[xTip,axisYZ];
preview.apexPosition=V(1,:);
preview.radialExtent=R;
preview.V=V; preview.F=F;
preview.meridianRadius=r;
preview.meridianX=x;
end
