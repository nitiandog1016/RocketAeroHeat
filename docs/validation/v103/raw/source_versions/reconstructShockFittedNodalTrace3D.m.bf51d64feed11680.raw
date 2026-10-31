function trace=reconstructShockFittedNodalTrace3D(Q,mesh,nodes,gas,nodeIDs,maximumArrayBytes)
%RECONSTRUCTSHOCKFITTEDNODALTRACE3D 在共享边界节点上拟合唯一的流体侧状态。
% 模块一：从节点关联四面体出发，只向内部邻接扩展；没有上游虚构样本。
% 模块二：距离加权的四未知仿射SVD拟合同时求节点值和XYZ梯度，保留秩与单元号。
% 模块三：以真实单元原始变量的凸均值为锚，显式共同斜率限制满足原物性域。
% 模块四：声明全部输入/输出和临时数组上界；不把插值校核当PDE或区域资格。
if nargin<6,maximumArrayBytes=128*1024^2;end
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:NodalTraceGas','节点重构当前仅常gamma。');
validateattributes(maximumArrayBytes,{'double'},{'scalar','finite','real','positive'});
assert(isequal(size(Q),[mesh.cellCount 5])&&isequal(size(nodes),[mesh.nodeCount 3])&& ...
    iscolumn(nodeIDs)&&~isempty(nodeIDs)&&all(isfinite(nodeIDs))&& ...
    all(nodeIDs==fix(nodeIDs)&nodeIDs>=1&nodeIDs<=mesh.nodeCount)&&numel(unique(nodeIDs))==numel(nodeIDs)&& ...
    all(ismember(nodeIDs,unique(mesh.faces(mesh.boundaryFaces,:)))), ...
    'ShockFit:NodalTraceShape','需要真实共享边界节点和对应单元五量。');
inputs=whos('Q','mesh','nodes','gas','nodeIDs');inputBytes=sum([inputs.bytes]);
temporaryBytes=8*(120*mesh.cellCount+20*mesh.faceCount+3*mesh.nodeCount+ ...
    numel(nodeIDs)*(100+4*mesh.cellCount));
budget=struct('retainedInputBytes',inputBytes,'declaredTemporaryAndOutputUpperBytes',temporaryBytes, ...
    'declaredTotalArrayUpperBytes',inputBytes+temporaryBytes,'maximumArrayBytes',maximumArrayBytes, ...
    'processPeakMemoryMeasured',false);
assert(budget.declaredTotalArrayUpperBytes<=maximumArrayBytes,'ShockFit:NodalTraceBudget', ...
    '共享节点重构的全部声明数组预算超限，未装配邻域。');
s=shockFittedState3D(Q,gas,'conservative');W=[s.density s.velocity s.pressure];
T=mesh.tetrahedra;centers=(nodes(T(:,1),:)+nodes(T(:,2),:)+nodes(T(:,3),:)+nodes(T(:,4),:))/4;
incident=sparse(T(:),repmat((1:mesh.cellCount)',4,1),1,mesh.nodeCount,mesh.cellCount);
inside=find(mesh.neighbour>0);a=mesh.owner(inside);b=mesh.neighbour(inside);
adjacency=sparse([a;b],[b;a],1,mesh.cellCount,mesh.cellCount);
count=numel(nodeIDs);stencils=cell(count,1);rings=zeros(count,1);singularValues=zeros(count,4);
anchor=zeros(count,5);increment=anchor;coefficient=zeros(4,5,count);lengthScale=zeros(count,1);
for j=1:count
    ids=find(incident(nodeIDs(j),:))';visited=false(mesh.cellCount,1);visited(ids)=true;front=ids;accepted=false;
    for ring=0:3
        ids=find(visited);
        if ring>=1&&numel(ids)>=6
            dx=centers(ids,:)-nodes(nodeIDs(j),:);distance=vecnorm(dx,2,2);h=max(distance);
            assert(h>0&&all(distance>0),'ShockFit:NodalTraceStencil','节点与单元中心重合。');
            center=mean(W(ids,:),1);
            if all(W(ids,:)==W(ids(1),:),'all'),center=W(ids(1),:);end
            A=[ones(numel(ids),1),dx/h]./distance;B=(W(ids,:)-center)./distance;
            [left,S,right]=svd(A,'econ');sigma=diag(S);
            if numel(sigma)==4&&sigma(4)>1e-10*sigma(1)
                fit=right*((left'*B)./sigma);
                anchor(j,:)=center;increment(j,:)=fit(1,:);coefficient(:,:,j)=fit;
                stencils{j}=ids;rings(j)=ring;singularValues(j,:)=sigma';lengthScale(j)=h;
                accepted=true;break;
            end
        end
        next=find(any(adjacency(front,:),1))';next=next(~visited(next));visited(next)=true;front=next;
    end
    assert(accepted,'ShockFit:NodalTraceRank','有限单侧节点邻域不能恢复四个仿射未知量。');
end
unlimited=anchor+increment;factor=ones(count,1);lower=gas.R*gas.temperatureRangeK(1);upper=gas.R*gas.temperatureRangeK(2);
base=[anchor(:,1),anchor(:,5),anchor(:,5)-lower*anchor(:,1),upper*anchor(:,1)-anchor(:,5)];
trial=[unlimited(:,1),unlimited(:,5),unlimited(:,5)-lower*unlimited(:,1),upper*unlimited(:,1)-unlimited(:,5)];
for j=1:count
    bad=trial(j,:)<0;bad(1:2)=trial(j,1:2)<=0;
    if any(bad),factor(j)=max(0,min([1 base(j,bad)./(base(j,bad)-trial(j,bad))]))*(1-64*eps);end
end
primitive=anchor+factor.*increment;temperature=primitive(:,5)./(gas.R*primitive(:,1));
assert(all(isfinite(primitive),'all')&&all(primitive(:,1)>0&primitive(:,5)>0)&& ...
    all(temperature>=gas.temperatureRangeK(1)&temperature<=gas.temperatureRangeK(2)), ...
    'ShockFit:NodalTraceAdmissibility','节点重构超出原物性域，不修改温度。');
state=shockFittedState3D([primitive(:,1:4) temperature],gas,'primitive');nodeQ=state.Q;
for j=1:count
    ids=stencils{j};
    if all(Q(ids,:)==Q(ids(1),:),'all'),nodeQ(j,:)=Q(ids(1),:);end
end
trace=struct('nodeIDs',nodeIDs,'nodePositionM',nodes(nodeIDs,:),'cellCentroidM',centers, ...
    'stencilCellIDs',{stencils},'stencilRings',rings,'stencilSingularValues',singularValues, ...
    'scaledAffineCoefficientsRP',coefficient,'fitLengthScaleM',lengthScale, ...
    'referencePrimitiveRP',anchor,'unlimitedPrimitiveRP',unlimited,'primitiveRP',primitive, ...
    'slopeFactor',factor,'limitedNodeCount',nnz(factor<1),'Q',nodeQ,'budget',budget, ...
    'primitiveConvention','[rho ux uy uz p]；仿射系数行分别为常数及归一化XYZ', ...
    'method','真实节点关联流体单元及有限邻接环/距离加权四未知SVD/共同物性域斜率系数', ...
    'cellConservedStateChanged',false,'externalCFDImported',false,'productionQualified',false);
end
