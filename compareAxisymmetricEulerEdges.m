function report = compareAxisymmetricEulerEdges(first,second)
%COMPAREAXISYMMETRICEULEREDGES 对比同外形、同工况、同算法的两个网格。
% 模块用途：检查近壁状态及沿程梯度对网格是否敏感；不当作真实误差。
% 只接受已收敛且远场合格的场。保持第一套壁面中心为公共采样位置，
% 第二套结果沿母线分段线性插值；不外推，也不越过折角或逆流段。
% 两网格一致并不足以证明物理解正确，因此不自动开放热流入口。
if ~first.available || ~second.available || ~first.converged || ~second.converged || ...
        ~first.boundary.clear || ~second.boundary.clear
    error('RocketAeroHeat:EulerCompare:Unconverged','两个场均须收敛且通过远场检查。');
end
for name={'gamma','R','Mach','speed','temperature','pressure','density','totalPressure','totalTemperature'}
    if ~isequal(first.freestream.(name{1}),second.freestream.(name{1}))
        error('RocketAeroHeat:EulerCompare:Flight','不能比较不同来流条件。');
    end
end
if ~strcmp(first.fluxScheme,second.fluxScheme)
    error('RocketAeroHeat:EulerCompare:Scheme','网格对比必须采用相同数值通量。');
end
orderOne=1;orderTwo=1;
if isfield(first,'spatialOrder'),orderOne=first.spatialOrder;end
if isfield(second,'spatialOrder'),orderTwo=second.spatialOrder;end
if orderOne~=orderTwo
    error('RocketAeroHeat:EulerCompare:Reconstruction','网格对比必须采用相同空间重构，不能把算法变化当网格误差。');
end
if orderOne==2
    for entry={{'boundaryReconstruction','oneSided'},{'reconstructionVariables','primitive'}}
        pair=entry{1};a=pair{2};b=pair{2};
        if isfield(first,pair{1}),a=first.(pair{1});end
        if isfield(second,pair{1}),b=second.(pair{1});end
        if ~strcmp(a,b)
            error('RocketAeroHeat:EulerCompare:Reconstruction','网格对比必须保持相同的边界重构和热力变量。');
        end
    end
end
if isequal(first.grid.X,second.grid.X) && isequal(first.grid.R,second.grid.R)
    error('RocketAeroHeat:EulerCompare:SameGrid','相同网格的续算不能当作网格加密对比。');
end
a=first.grid.contour;b=second.grid.contour;L=first.grid.length;
x=unique([a(:,1);b(:,1)]);
if a(1,1)~=b(1,1) || a(end,1)~=b(end,1) || ...
        max(abs(interp1(a(:,1),a(:,2),x)-interp1(b(:,1),b(:,2),x)))>1e-10*L
    error('RocketAeroHeat:EulerCompare:Geometry','两套网格的分片线性几何必须一致。');
end
one=extractAxisymmetricEulerEdge(first);two=extractAxisymmetricEulerEdge(second);
report.s=one.s;report.x=one.x;report.radius=one.radius;
report.method='相同几何/工况/通量的双网格近壁对比；公共弧长位置分段线性插值';
report.scope='只量化网格敏感性；不提供 Richardson 外推、真实误差或精度认证';
report.available=true;report.gridAccuracyVerified=false;report.edgeStateAvailable=false;
report.firstGrid=[first.grid.nx,first.grid.nr];report.secondGrid=[second.grid.nx,second.grid.nr];
report.axialNodesChanged=~isequal(first.grid.X(:,1),second.grid.X(:,1));
report.outerRadiusChanged=first.grid.outerRadius~=second.grid.outerRadius;
oneTopology='upstreamInlet';twoTopology='upstreamInlet';
if isfield(first.grid,'boundaryTopology'),oneTopology=first.grid.boundaryTopology;end
if isfield(second.grid,'boundaryTopology'),twoTopology=second.grid.boundaryTopology;end
report.boundaryTopologyChanged=~strcmp(oneTopology,twoTopology);
report.status='双网格差值不是独立精度证明；必须同时检查轴向、径向和几何分辨率';
if report.boundaryTopologyChanged
    report.status=[report.status '；本次还改变了远场/轴线边界拓扑，并非单纯单元加密'];
end
names={'pressure','temperature','density','signedSpeed','totalPressure','velocityGradient','pressureGradient'};
fs=first.freestream;
scales=[fs.pressure,fs.temperature,fs.density,fs.speed,fs.totalPressure,fs.speed/L,fs.pressure/L];
mask=~one.flags.feature & ~one.flags.reverseFlow;
for k=1:numel(names)
    name=names{k};v1=one.(name);
    v2=interpolateSegments(two.s,two.(name), ...
        ~two.flags.feature & ~two.flags.reverseFlow,one.s);
    valid=mask & isfinite(v1) & isfinite(v2);
    if nnz(valid)<3
        error('RocketAeroHeat:EulerCompare:Overlap','公共光滑区采样不足，不能比较 %s。',name);
    end
    delta=v2-v1;delta(~valid)=NaN;w=one.segmentLength(valid);w=w/sum(w);
    metric=struct('first',v1,'secondInterpolated',v2,'difference',delta,'valid',valid, ...
        'referenceScale',scales(k),'maximumAbsoluteDifference',max(abs(delta(valid))), ...
        'rmsDifference',sqrt(sum(w.*delta(valid).^2)), ...
        'maximumNormalizedDifference',max(abs(delta(valid)))/scales(k));
    report.(name)=metric;
end
end

function values=interpolateSegments(s,y,valid,target)
% 模块作用：逐段插值，保留折角/逆流/无导数区造成的缺口。
valid=valid & isfinite(y);values=nan(size(target));
starts=find(valid & [true;~valid(1:end-1)]);ends=find(valid & [~valid(2:end);true]);
for j=1:numel(starts)
    ids=starts(j):ends(j);
    if numel(ids)<2,continue;end
    selected=target>=s(ids(1)) & target<=s(ids(end));
    values(selected)=interp1(s(ids),y(ids),target(selected),'linear');
end
end
