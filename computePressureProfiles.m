function profiles = computePressureProfiles(mesh)
%COMPUTEPRESSUREPROFILES 从现有面元压力提取柱段周向、轴向诊断曲线。
% 不插值生成新的压力解；每个点是对应角度或轴向区间的面积加权平均。
% 仅选取法向 X 分量较小的近圆柱侧壁，避免鼻部高压混入柱段曲线。
% 输出包含所用截面位置、角度窗口和每个采样点的面元数量。

if ~isfield(mesh,'pressure') || ~isfield(mesh,'localGeometry')
    error('RocketAeroHeat:Profiles:Input','请先计算局部几何与压力。');
end
g=mesh.localGeometry; p=mesh.pressure;
nf=size(mesh.F,1);
if numel(p.Cp)~=nf || numel(g.theta)~=nf
    error('RocketAeroHeat:Profiles:Input','压力和表面坐标与网格不匹配。');
end

%% 1. 识别近圆柱段；不足时明确报告，避免生成误导性曲线
mask=g.isSidewall(:) & abs(g.normal(:,1))<0.02 & ...
    isfinite(g.theta(:)) & isfinite(g.area(:)) & g.area(:)>0;
profiles=struct('available',false,'reason','','cylinderMask',mask);
if nnz(mask)<50
    profiles.reason='近圆柱面元不足，无法生成柱段压力曲线';
    return;
end
x=g.center(:,1); theta=mod(rad2deg(g.theta(:)),360);
cp=p.Cp(:); area=g.area(:);
xValues=x(mask); xMin=min(xValues); xMax=max(xValues);
if xMax-xMin<=0
    profiles.reason='柱段轴向跨度为零';
    return;
end

%% 2. 中间截面的周向曲线：在有限宽度轴向带内分角度箱求均值
stationX=(xMin+xMax)/2;
width=max(0.015*(max(mesh.V(:,1))-min(mesh.V(:,1))),0.02*(xMax-xMin));
station=mask & abs(x-stationX)<=width/2;
if nnz(station)<24
    width=2*width; station=mask & abs(x-stationX)<=width/2;
end
if nnz(station)<24
    profiles.reason='选定截面的周向采样不足';
    return;
end
angleEdges=0:10:360;
[angularCp,angularCount]=binWeighted(theta(station),cp(station),area(station),angleEdges);

%% 3. 沿轴向统计迎风线与背风线：窗口中心随当前来流改变
flow=p.freestream.unitFlowDirection(:)';
if norm(flow(2:3))<1e-10
    windwardTheta=0; % 零攻角轴对称，选 +Y 方向作参考。
else
    windwardTheta=mod(atan2d(-flow(3),-flow(2)),360);
end
leewardTheta=mod(windwardTheta+180,360);
windowDeg=15;
windward=mask & angleDistance(theta,windwardTheta)<=windowDeg;
leeward=mask & angleDistance(theta,leewardTheta)<=windowDeg;
xEdges=linspace(xMin,xMax,81);
[windCp,windCount]=binWeighted(x(windward),cp(windward),area(windward),xEdges);
[leeCp,leeCount]=binWeighted(x(leeward),cp(leeward),area(leeward),xEdges);

%% 4. 返回原始物理坐标、采样数量和用于 GUI 色标的柱段范围
profiles.available=true; profiles.reason='';
profiles.stationX=stationX; profiles.stationWidth=width;
profiles.stationFaceCount=nnz(station);
profiles.angleDeg=(angleEdges(1:end-1)+angleEdges(2:end))'/2;
profiles.angularCp=angularCp; profiles.angularCount=angularCount;
profiles.axialX=(xEdges(1:end-1)+xEdges(2:end))'/2;
profiles.windwardCp=windCp; profiles.windwardCount=windCount;
profiles.leewardCp=leeCp; profiles.leewardCount=leeCount;
profiles.windwardThetaDeg=windwardTheta;
profiles.leewardThetaDeg=leewardTheta;
profiles.angularWindowDeg=windowDeg;
profiles.cylinderCpRange=[min(cp(mask)) max(cp(mask))];
end

function distance=angleDistance(theta,center)
% 角度差折回 [-180°,180°]，保证 0°/360° 接缝连续。
distance=abs(mod(theta-center+180,360)-180);
end

function [average,count]=binWeighted(location,value,weight,edges)
% 每个区间按面元面积平均；空区间设为 NaN，不伪造连续数据。
[~,~,bin]=histcounts(location,edges);
valid=bin>0 & isfinite(value) & isfinite(weight) & weight>0;
n=numel(edges)-1;
count=accumarray(bin(valid),1,[n 1],@sum,0);
weighted=accumarray(bin(valid),value(valid).*weight(valid),[n 1],@sum,0);
totalWeight=accumarray(bin(valid),weight(valid),[n 1],@sum,0);
average=nan(n,1); occupied=totalWeight>0;
average(occupied)=weighted(occupied)./totalWeight(occupied);
end
