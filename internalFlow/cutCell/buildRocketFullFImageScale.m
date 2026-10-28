function [columnScale,columnNorm,report]=buildRocketFullFImageScale(graph,factor,options)
%BUILDROCKETFULLFIMAGESCALE 流式计算完整A_FF*M_F每个像列的二范数，不形成全局乘积。
% 模块一：按原列建立入边索引；全部原F五量列按固定局部组恰好覆盖一次。
% 模块二：每组收集所有真实受影响行，乘已认证的实际局部逆，包含全部跨组连接。
% 模块三：严格正范数取倒数作右侧坐标尺度，不裁原系数/方程或改变原残差范数。
if nargin<3,options=struct;end
defaults=struct('maximumWorkspaceBytes',256*1024^2);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'FImageScale:Options','未知完整像列归一选项。');
if ~isfield(options,'maximumWorkspaceBytes'),options.maximumWorkspaceBytes=defaults.maximumWorkspaceBytes;end
assert(strcmp(graph.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    strcmp(factor.kind,'fixed-within-original-patch-five-field-F-inverse')&&graph.cellCount==factor.cellCount,'FImageScale:Input','需要原完整F图及其固定局部逆。');
n=graph.cellCount;nb=numel(graph.column);K=numel(factor.buckets);counts=accumarray(double(graph.column),1,[n 1]);
% 此处是额外算法数组预算，不含调用方的原图/因子，也不是进程峰值实测。
bound=32*nb+128*n+16*25*K^2*max(counts)+65536;
assert(isfinite(options.maximumWorkspaceBytes)&&options.maximumWorkspaceBytes>0&&bound<=options.maximumWorkspaceBytes,'FImageScale:Budget','完整像列归一额外数组超预算，拒绝而不删列。');
[~,entries]=sort(graph.column);pointer=[1;1+cumsum(counts)];owners=repelem((1:n)',diff(double(graph.rowPointer)));
columnNorm=zeros(n,5);coverage=zeros(n,1);used=0;groups=0;maximumPanelBytes=0;maximumRows=0;timer=tic;
for k=1:K
    bucket=factor.buckets{k};
    for page=1:numel(bucket.groupIDs)
        physical=double(bucket.members(:,page));algebraic=double(graph.inverseOrder(physical));selected=[];localColumn=[];
        for j=1:k
            ids=double(entries(pointer(algebraic(j)):pointer(algebraic(j)+1)-1));
            selected=[selected;ids(:)];localColumn=[localColumn;repmat(j,numel(ids),1)]; %#ok<AGROW>
        end
        [rows,~,localRow]=unique(owners(selected));panel=zeros(5*numel(rows),5*k);
        for a=1:5,for b=1:5
            index=sub2ind(size(panel),5*(localRow-1)+a,5*(localColumn-1)+b);panel(index)=reshape(graph.values(a,b,selected),[],1);
        end,end
        image=panel*bucket.inverse(:,:,page);norms=vecnorm(image,2,1);
        assert(all(isfinite(norms)&norms>0),'FImageScale:ZeroColumn','完整预条件像列非有限或为零，不补小数或截断。');
        assert(all(coverage(physical)==0),'FImageScale:Coverage','原F列重复覆盖。');
        columnNorm(physical,:)=reshape(norms,5,k)';coverage(physical)=coverage(physical)+1;
        used=used+numel(selected);groups=groups+1;maximumRows=max(maximumRows,numel(rows));
        p=whos('panel','image');maximumPanelBytes=max(maximumPanelBytes,sum([p.bytes]));
    end
end
assert(all(coverage==1)&&used==nb,'FImageScale:Coverage','完整F五量列或原块未完全且唯一覆盖。');
columnScale=1./columnNorm;assert(all(isfinite(columnScale)&columnScale>0,'all'),'FImageScale:Scale','像列倒数非有限或非正。');
report=struct('method','inverse-complete-A_FF-M_F-column-two-norm','groups',groups,'originalBlocksUsedExactlyOnce',used, ...
    'allOriginalUnknownColumnsCoveredExactlyOnce',true,'columnNormRange',[min(columnNorm,[],'all') max(columnNorm,[],'all')], ...
    'columnNormSpan',max(columnNorm,[],'all')/min(columnNorm,[],'all'),'columnScaleRange',[min(columnScale,[],'all') max(columnScale,[],'all')], ...
    'maximumSupportedOriginalRows',maximumRows,'maximumDensePanelAndImageArrayBytes',maximumPanelBytes, ...
    'additionalWorkspaceArrayBoundBytes',bound,'fullGlobalPreconditionedMatrixAllocated',false,'processPeakMemoryClaimed',false, ...
    'fixedPositiveRightScalingOnly',true,'originalFResidualNormAndRightHandSideUnchanged',true,'elapsedSeconds',toc(timer), ...
    'newLocalInverseConstructions',0,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
