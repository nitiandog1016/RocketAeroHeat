function [profile,stations] = computeCircumferentialPressure(mesh,targetX,stations)
%COMPUTECIRCUMFERENTIALPRESSURE 任意轴向位置的周向压力诊断。
% 按现有网格的相邻截面环带统计；没有截面拓扑时按轴向小区间统计。
% 每个角度点是该位置面元的面积加权平均，空角度区间为 NaN。
% 输入 targetX 为模型坐标 X（m）；输出 faceIDs 可用于三维标记。
% stations 可复用，拖动滑块时无需每次重新建立环带列表。
if ~all(isfield(mesh,{'V','F','pressure','localGeometry','surfaceCoord'}))
    error('RocketAeroHeat:Circumference:Input','需要最终网格、表面坐标和压力。');
end
nf=size(mesh.F,1); g=mesh.localGeometry; sc=mesh.surfaceCoord;
if numel(sc.faceTheta)~=nf || numel(mesh.pressure.Cp)~=nf || ...
        numel(g.area)~=nf || numel(g.center(:,1))~=nf
    error('RocketAeroHeat:Circumference:Size','面元压力或表面坐标与网格不匹配。');
end
validateattributes(targetX,{'numeric'},{'real','finite','scalar'});
if nargin<3 || isempty(stations)
    stations=makeStations(mesh);
end
if isempty(stations.x)
    error('RocketAeroHeat:Circumference:NoStation','模型缺少可用的轴向截面。');
end

%% 1. 选距离目标 X 最近的实际环带，显示该环带的真实轴向范围
[~,index]=min(abs(stations.x-targetX));
faceIDs=stations.faceIDs{index};
theta=mod(rad2deg(sc.faceTheta(faceIDs)),360);
area=g.area(faceIDs); cp=mesh.pressure.Cp(faceIDs);
pressure=mesh.pressure.staticPressure(faceIDs);
valid=isfinite(theta) & isfinite(area) & area>0 & ...
    isfinite(cp) & isfinite(pressure);
theta=theta(valid); area=area(valid);
cp=cp(valid); pressure=pressure(valid);
faceIDs=faceIDs(valid);

%% 2. 将周向角分成 36 个 10° 区间；空箱保持 NaN，不插值补点
edges=0:10:360; [~,~,bin]=histcounts(theta,edges); keep=bin>0;
bin=bin(keep); area=area(keep); cp=cp(keep); pressure=pressure(keep);
count=accumarray(bin,1,[36 1],@sum,0);
weight=accumarray(bin,area,[36 1],@sum,0);
cpSum=accumarray(bin,cp.*area,[36 1],@sum,0);
pSum=accumarray(bin,pressure.*area,[36 1],@sum,0);
cpMean=nan(36,1); pMean=nan(36,1);
occupied=weight>0;
cpMean(occupied)=cpSum(occupied)./weight(occupied);
pMean(occupied)=pSum(occupied)./weight(occupied);

%% 3. 保存采样含义，便于解释尖鼻与封盖处较少的周向数据
profile=struct;
profile.requestedX=targetX; profile.stationX=stations.x(index);
profile.xRange=stations.xRange(index,:);
profile.stationIndex=index; profile.stationCount=numel(stations.x);
profile.source=stations.source{index};
profile.faceIDs=faceIDs; profile.faceCount=numel(faceIDs);
profile.angleDeg=(5:10:355)';
profile.Cp=cpMean; profile.staticPressurePa=pMean;
profile.count=count; profile.coverage=nnz(occupied)/36;
profile.available=nnz(occupied)>=3;
if profile.available
    profile.message='';
else
    profile.message='此处周向面元不足，曲线仅显示已有角度或留空';
end
end

function stations=makeStations(mesh)
% 优先复用网格原生 faceBandID；平封盖单独成为端部站位。
nf=size(mesh.F,1); center=mesh.localGeometry.center(:,1);
area=mesh.localGeometry.area(:); usable=isfinite(center) & isfinite(area) & area>0;
groups={}; source={};
assigned=false(nf,1);
if isfield(mesh.surfaceCoord,'faceBandID') && numel(mesh.surfaceCoord.faceBandID)==nf
    band=mesh.surfaceCoord.faceBandID(:);
    ids=unique(band(usable & isfinite(band) & band>0));
    for j=1:numel(ids)
        members=find(usable & band==ids(j));
        if isempty(members), continue; end
        groups{end+1,1}=members; %#ok<AGROW>
        source{end+1,1}='原始轴向环带'; %#ok<AGROW>
        assigned(members)=true;
    end
end
% 两端平封盖有相同的 bandID 时仍分开，避免首尾混为一个截面。
for key={'isAftCap','isNoseCap'}
    if isfield(mesh.localGeometry,key{1})
        capFlag=mesh.localGeometry.(key{1});
        members=find(usable & capFlag(:));
        members=members(~assigned(members));
        if ~isempty(members)
            groups{end+1,1}=members; %#ok<AGROW>
            source{end+1,1}='端部封盖面元'; %#ok<AGROW>
            assigned(members)=true;
        end
    end
end
% 普通导入网格缺少环带信息时，按轴向位置分成最多 120 个取样区。
remaining=find(usable & ~assigned);
if ~isempty(remaining)
    x=center(remaining); width=max(x)-min(x);
    if width==0
        bins=ones(numel(remaining),1); nbin=1;
    else
        nbin=min(120,max(1,round(numel(remaining)/30)));
        bins=min(nbin,1+floor((x-min(x))/width*nbin));
    end
    for j=1:nbin
        members=remaining(bins==j);
        if isempty(members), continue; end
        groups{end+1,1}=members; %#ok<AGROW>
        source{end+1,1}='轴向面元分箱'; %#ok<AGROW>
    end
end
nx=numel(groups); x=zeros(nx,1); limits=zeros(nx,2);
for j=1:nx
    members=groups{j}; localX=center(members); weights=area(members);
    x(j)=sum(localX.*weights)/sum(weights);
    limits(j,:)=[min(localX) max(localX)];
end
[x,order]=sort(x); stations=struct;
stations.x=x; stations.xRange=limits(order,:);
stations.faceIDs=groups(order); stations.source=source(order);
end
