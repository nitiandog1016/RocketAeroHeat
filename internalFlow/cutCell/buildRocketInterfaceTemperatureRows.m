function rows=buildRocketInterfaceTemperatureRows(Q,cut,gas,scale,faceIDs)
%BUILDROCKETINTERFACETEMPERATUREROWS 构造真实内部Riemann采样温度对缩放守恒方向的解析行。
% 模块一：显式选择当前内部面，调用现有精确分支采样导数，保留原温域和波系护栏。
% 模块二：每面只连接两侧原CV的十个未知量，原CV/五量排列保持，不补物理边界条件。
% 模块三：这些行只用于研究迭代方向约束；不能保证非线性候选温度或授予区域资格。
n=size(Q,1);validateattributes(faceIDs,{'numeric'},{'vector','integer','positive','finite','nonempty'});
faceIDs=faceIDs(:);assert(numel(unique(faceIDs))==numel(faceIDs), ...
    'TemperatureRows:Faces','显式内部面不得重复，不隐式删约束。');
m=numel(faceIDs);
% 模块四：按16面调用原导数内核，全部元数据和十项稀疏导数保留。
% 64MiB约束为额外构造数组预算，不含输入Q/几何；不按面数分配整箭稠密行。
assert(8*(5*n+1)+2048*m+65536<=64*1024^2, ...
    'TemperatureRows:Budget','全部约束行的额外构造数组超过64MiB，不截断约束。');
assert(isequal(size(scale),[1 5])&&all(isfinite(scale)&scale>0),'TemperatureRows:Scale','需要真实五量正缩放。');
f=cut.openFaces;assert(all(faceIDs<=numel(f.area))&&all(f.area(faceIDs)>0&f.cellMinus(faceIDs)>0&f.cellPlus(faceIDs)>0), ...
    'TemperatureRows:Faces','温度约束只接受当前实际内部面。');
owner=double(f.cellMinus(faceIDs));neighbour=double(f.cellPlus(faceIDs));normal=f.normal(faceIDs,:);
QL=Q(owner,:);QR=Q(neighbour,:);
left=zeros(m,5);right=left;basePrimitive=zeros(m,5);sampleRegion=zeros(m,1);
for first=1:16:m
    ids=first:min(first+15,m);count=numel(ids);
    [~,~,audit,sample]=cutCellExactEulerBranchJacobian(QL(ids,:),QR(ids,:),normal(ids,:),gas);
    left(ids,:)=reshape(sample.left(5,:,:),5,count)'.*scale;
    right(ids,:)=reshape(sample.right(5,:,:),5,count)'.*scale;
    basePrimitive(ids,:)=sample.basePrimitive;sampleRegion(ids)=audit.sampleRegion;
end
columns=[owner+(0:4)*n,neighbour+(0:4)*n];values=[left right];
rowIndex=repmat((1:m)',1,10);H=sparse(rowIndex(:),columns(:),values(:),m,5*n);
norms=sqrt(full(sum(H.^2,2)));assert(all(isfinite(norms)&norms>0),'TemperatureRows:Derivative','采样温度导数必须有限且非零。');
rows=struct('kind','current-interior-selected-Riemann-sample-temperature-tangent-rows','cellCount',n, ...
    'faceIDs',faceIDs,'ownerCV',owner,'neighbourCV',neighbour,'normal',normal,'QL',QL,'QR',QR, ...
    'scale',scale,'gas',gas,'H',H,'leftRowsK',left,'rightRowsK',right, ...
    'baseSamplePrimitive',basePrimitive,'baseTemperatureK',basePrimitive(:,5), ...
    'remainingUpperTemperatureK',gas.temperatureRange(2)-basePrimitive(:,5), ...
    'sampleRegion',sampleRegion,'rowTwoNormK',norms,'sourceStateChanged',false, ...
    'nonlinearTemperatureBoundGuaranteed',false,'physicalBoundaryConditionAdded',false, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
