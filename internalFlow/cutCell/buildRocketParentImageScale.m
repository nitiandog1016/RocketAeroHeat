function [columnScale,columnNorm,report]=buildRocketParentImageScale(packed,factor,options)
%BUILDROCKETPARENTIMAGESCALE 流式求完整父C乘原五量局部映射的全部像列范数。
% 模块一：完整父矩阵的代数行列和原CV排列均保留，不使用C/F消元或旧Schur。
% 模块二：每组收集所有受影响原行，乘缓存实际局部逆；每个原块恰读取一次。
% 模块三：仅取严格正完整列范数的倒数作右侧尺度，原方程和残差范数不变。
if nargin<3,options=struct;end
defaults=struct('maximumWorkspaceBytes',256*1024^2);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'ParentImageScale:Options','未知完整父矩阵像列选项。');
if ~isfield(options,'maximumWorkspaceBytes'),options.maximumWorkspaceBytes=defaults.maximumWorkspaceBytes;end
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    strcmp(factor.kind,'fixed-five-conservation-strong-connection-local-block-inverse')&&packed.cellCount==factor.cellCount,'ParentImageScale:Input','需要完整父矩阵及全部原五量局部映射。');
n=packed.cellCount;nb=numel(packed.column);K=numel(factor.buckets);counts=accumarray(double(packed.column),1,[n 1]);
bound=32*nb+128*n+16*25*K^2*max(counts)+65536;
assert(isfinite(options.maximumWorkspaceBytes)&&options.maximumWorkspaceBytes>0&&bound<=options.maximumWorkspaceBytes,'ParentImageScale:Budget','全部像列额外数组超预算，不删列。');
[~,entries]=sort(packed.column);pointer=[1;1+cumsum(counts)];owners=repelem((1:n)',diff(double(packed.rowPointer)));
columnNorm=zeros(n,5);coverage=zeros(n,1);used=0;groups=0;maximumRows=0;maximumPanelBytes=0;timer=tic;
for k=1:K
    bucket=factor.buckets{k};
    for page=1:numel(bucket.patchIDs)
        physical=double(bucket.members(:,page));algebraic=double(packed.inverseOrder(physical));selected=[];localColumn=[];
        for j=1:k
            ids=double(entries(pointer(algebraic(j)):pointer(algebraic(j)+1)-1));
            selected=[selected;ids(:)];localColumn=[localColumn;repmat(j,numel(ids),1)]; %#ok<AGROW>
        end
        [rows,~,localRow]=unique(owners(selected));panel=zeros(5*numel(rows),5*k);
        for a=1:5,for b=1:5
            index=sub2ind(size(panel),5*(localRow-1)+a,5*(localColumn-1)+b);panel(index)=reshape(packed.values(a,b,selected),[],1);
        end,end
        image=panel*bucket.inverse(:,:,page);norms=vecnorm(image,2,1);
        assert(all(isfinite(norms)&norms>0),'ParentImageScale:ZeroColumn','完整父矩阵预条件像列非有限或为零，不裁剪或补小数。');
        assert(all(coverage(physical)==0),'ParentImageScale:Coverage','原CV像列重复覆盖。');
        columnNorm(physical,:)=reshape(norms,5,k)';coverage(physical)=1;used=used+numel(selected);groups=groups+1;
        p=whos('panel','image');maximumPanelBytes=max(maximumPanelBytes,sum([p.bytes]));maximumRows=max(maximumRows,numel(rows));
    end
end
assert(all(coverage==1)&&used==nb,'ParentImageScale:Coverage','全部原CV五量像列或原块覆盖不全。');
columnScale=1./columnNorm;assert(all(isfinite(columnScale)&columnScale>0,'all'),'ParentImageScale:Scale','完整像列倒数非有限或非正。');
report=struct('method','inverse-complete-parent-C-local-five-field-map-column-two-norm','groups',groups,'originalBlocksUsedExactlyOnce',used, ...
    'allOriginalCVs',n,'allOriginalUnknownColumnsCoveredExactlyOnce',true,'columnNormRange',[min(columnNorm,[],'all') max(columnNorm,[],'all')], ...
    'columnNormSpan',max(columnNorm,[],'all')/min(columnNorm,[],'all'),'columnScaleRange',[min(columnScale,[],'all') max(columnScale,[],'all')], ...
    'maximumSupportedOriginalRows',maximumRows,'maximumDensePanelAndImageArrayBytes',maximumPanelBytes, ...
    'additionalWorkspaceArrayBoundBytes',bound,'fullGlobalPreconditionedMatrixAllocated',false,'processPeakMemoryClaimed',false, ...
    'oldSchurUsed',false,'fixedCRepresentativesUsed',false,'originalParentResidualNormUnchanged',true,'elapsedSeconds',toc(timer), ...
    'newLocalInverseConstructions',0,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
