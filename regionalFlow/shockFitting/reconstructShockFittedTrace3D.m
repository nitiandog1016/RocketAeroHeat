function trace=reconstructShockFittedTrace3D(Q,mesh,nodes,gas,faceIDs,options)
%RECONSTRUCTSHOCKFITTEDTRACE3D 从流体侧单元中心重构真实边界面上的下游状态。
% 模块一：只用内部邻接单元，按距离加权最小二乘恢复[rho,u,v,w,p]三维梯度。
% 模块二：必要时扩大有限环数，检查三维秩；缺失方向不能以零梯度冒充重构。
% 模块三：可显式启用共同斜率系数，满足密度/压力/原温域的线性不等式。
% 模块四：保留未限制试探值、系数、全部单元号与梯度，不改写单元五量或放宽温域。
% 模块五：推进前声明输入/输出/临时数组预算；它不是MATLAB进程峰值。
if nargin<6,options=struct;end
defaults=struct('initialRings',2,'maximumRings',3,'minimumNeighbors',6, ...
    'slopeLimiter','admissibility','maximumArrayBytes',128*1024^2);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'ShockFit:TraceOptions','存在未知的边界重构选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.initialRings,{'double'},{'scalar','real','finite','integer','>=',1,'<=',4});
validateattributes(options.maximumRings,{'double'},{'scalar','real','finite','integer','>=',options.initialRings,'<=',4});
validateattributes(options.minimumNeighbors,{'double'},{'scalar','real','finite','integer','>=',3});
validateattributes(options.maximumArrayBytes,{'double'},{'scalar','real','finite','positive'});
assert(ischar(options.slopeLimiter)&&ismember(options.slopeLimiter,{'none','admissibility'}), ...
    'ShockFit:TraceOptions','斜率限制必须显式为none或admissibility。');
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:TraceGas', ...
    '本三维界面重构仅常gamma，NASA7没有实现。');
assert(isequal(size(Q),[mesh.cellCount 5])&&isequal(size(nodes),[mesh.nodeCount 3]), ...
    'ShockFit:TraceShape','五量与当前共享节点形状不符。');
assert(iscolumn(faceIDs)&&~isempty(faceIDs)&&all(isfinite(faceIDs))&& ...
    all(faceIDs==fix(faceIDs)&faceIDs>=1&faceIDs<=mesh.faceCount)&& ...
    numel(unique(faceIDs))==numel(faceIDs)&&all(mesh.neighbour(faceIDs)==0), ...
    'ShockFit:TraceFaces','须提供不重复的真实边界面号，不能跨激波取外侧单元。');
input=whos('Q','mesh','nodes','gas','faceIDs','options');inputBytes=sum([input.bytes]);
maximumStencil=min(mesh.cellCount-1,2*(3^options.maximumRings-1));
temporaryBytes=8*(100*mesh.cellCount+20*mesh.faceCount+3*mesh.nodeCount+ ...
    numel(faceIDs)*(60+4*maximumStencil));
budget=struct('retainedInputBytes',inputBytes,'declaredTemporaryAndOutputUpperBytes',temporaryBytes, ...
    'declaredTotalArrayUpperBytes',inputBytes+temporaryBytes,'maximumArrayBytes',options.maximumArrayBytes, ...
    'maximumStencilCellsPerOwner',maximumStencil,'processPeakMemoryMeasured',false);
assert(budget.declaredTotalArrayUpperBytes<=options.maximumArrayBytes,'ShockFit:TraceMemoryBudget', ...
    '边界重构的声明数组预算超限，尚未装配邻接或梯度。');
state=shockFittedState3D(Q,gas,'conservative');W=[state.density state.velocity state.pressure];
T=mesh.tetrahedra;cellCenters=(nodes(T(:,1),:)+nodes(T(:,2),:)+nodes(T(:,3),:)+nodes(T(:,4),:))/4;
F=mesh.faces(faceIDs,:);faceCenters=(nodes(F(:,1),:)+nodes(F(:,2),:)+nodes(F(:,3),:))/3;
inside=find(mesh.neighbour>0);a=mesh.owner(inside);b=mesh.neighbour(inside);
adjacency=sparse([a;b],[b;a],1,mesh.cellCount,mesh.cellCount);
owners=mesh.owner(faceIDs);[uniqueOwners,~,ownerMap]=unique(owners);count=numel(uniqueOwners);
gradient=zeros(3,5,count);stencil=cell(count,1);rings=zeros(count,1);singularValues=zeros(count,3);
for j=1:count
    owner=uniqueOwners(j);visited=false(mesh.cellCount,1);visited(owner)=true;front=owner;accepted=false;
    for ring=1:options.maximumRings
        next=find(any(adjacency(front,:),1))';next=next(~visited(next));visited(next)=true;front=next;
        ids=find(visited);ids(ids==owner)=[];
        if ring<options.initialRings||numel(ids)<options.minimumNeighbors,continue;end
        dx=cellCenters(ids,:)-cellCenters(owner,:);distance=vecnorm(dx,2,2);
        assert(all(distance>0),'ShockFit:TraceStencil','不同单元中心重合，距离权重无定义。');
        A=dx./distance;B=(W(ids,:)-W(owner,:))./distance;
        [left,S,right]=svd(A,'econ');s=diag(S);
        if numel(s)<3||s(3)<=1e-10*s(1),continue;end
        gradient(:,:,j)=right*((left'*B)./s);stencil{j}=ids;rings(j)=ring;singularValues(j,:)=s';
        accepted=true;break;
    end
    assert(accepted,'ShockFit:TraceStencilRank','有限单侧邻域不能恢复三维梯度，不补缺失方向。');
end
increment=zeros(numel(faceIDs),5);
for j=1:numel(faceIDs)
    increment(j,:)=(faceCenters(j,:)-cellCenters(owners(j),:))*gradient(:,:,ownerMap(j));
end
center=W(owners,:);unlimited=center+increment;factor=ones(numel(faceIDs),1);
if strcmp(options.slopeLimiter,'admissibility')
    % T=p/(rho*R)的原上下界，在rho>0时均是rho/p的线性不等式。
    lower=gas.R*gas.temperatureRangeK(1);upper=gas.R*gas.temperatureRangeK(2);
    base=[center(:,1),center(:,5),center(:,5)-lower*center(:,1),upper*center(:,1)-center(:,5)];
    trial=[unlimited(:,1),unlimited(:,5),unlimited(:,5)-lower*unlimited(:,1),upper*unlimited(:,1)-unlimited(:,5)];
    assert(all(base(:,1:2)>0,'all')&&all(base(:,3:4)>=-64*eps*max(center(:,5)),'all'), ...
        'ShockFit:TraceAdmissibility','原单元状态不满足重构的原物性域。');
    for j=1:numel(faceIDs)
        bad=trial(j,:)<0;bad(1:2)=trial(j,1:2)<=0;
        if any(bad)
            ratio=base(j,bad)./(base(j,bad)-trial(j,bad));
            factor(j)=max(0,min([1 ratio]))*(1-64*eps);
        end
    end
end
reconstructed=center+factor.*increment;
temperature=reconstructed(:,5)./(gas.R*reconstructed(:,1));
assert(all(isfinite(reconstructed),'all')&&all(reconstructed(:,1)>0&reconstructed(:,5)>0)&& ...
    all(temperature>=gas.temperatureRangeK(1)&temperature<=gas.temperatureRangeK(2)), ...
    'ShockFit:TraceAdmissibility','面上重构超出原物性域；未修改温度或单元五量。');
faceState=shockFittedState3D([reconstructed(:,1:4) temperature],gas,'primitive');
faceQ=faceState.Q;unchanged=all(increment==0,2)|factor==0;
faceQ(unchanged,:)=Q(owners(unchanged),:);
trace=struct('faceIDs',faceIDs,'ownerCells',owners,'cellCentroidM',cellCenters, ...
    'faceCentroidM',faceCenters,'uniqueOwnerCells',uniqueOwners,'ownerGradientIndex',ownerMap, ...
    'gradientRP',gradient,'stencilCellIDs',{stencil},'stencilRings',rings, ...
    'stencilSingularValues',singularValues,'cellPrimitiveRP',center, ...
    'unlimitedPrimitiveRP',unlimited,'primitiveRP',reconstructed,'slopeFactor',factor, ...
    'limitedFaceCount',nnz(factor<1),'Q',faceQ,'budget',budget, ...
    'primitiveConvention','[rho ux uy uz p]；梯度第1维为XYZ，温度由原状态方程恢复', ...
    'method','流体侧有限环距离加权SVD最小二乘/面重心线性外推/显式共同斜率系数', ...
    'slopeLimiter',options.slopeLimiter,'cellConservedStateChanged',false, ...
    'externalCFDImported',false,'shockStatePrescribed',false,'productionQualified',false);
end
