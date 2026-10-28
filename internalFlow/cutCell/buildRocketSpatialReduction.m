function reduction=buildRocketSpatialReduction(packed,partition,coarseOriginalCVs,options)
%BUILDROCKETSPATIALREDUCTION 构造空间五量近似Schur S=D-E*M_F*B，保留所有填充。
% 模块一：F/C各保留整个原CV的五量，M_F仅为认证的原组内F逆。
% 模块二：所有原F-C/C-F/C-C块及局部消元生成的位置完整保留，不补对角。
% 模块三：本构造不求粗逆，不声称S是完整A_FF精确消元所得的Schur。
if nargin<4,options=struct;end
defaults=struct('maximumReductionBytes',768*1024^2,'maximumFFFactorBytes',384*1024^2,'maximumConstructionBytes',3072*1024^2, ...
    'maximumSchurBlocks',1000000,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'SpatialReduction:Options','未知空间消元选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
if isempty(coarseOriginalCVs),split=buildRocketSpatialSplit(packed,partition);else,split=buildRocketSpatialSplit(packed,partition,coarseOriginalCVs);end
n=split.cellCount;nc=split.coarseCells;nf=split.fineCells;nb=numel(packed.column);fCounts=double(split.groupFCounts);
factorBound=8*25*sum(fCounts.^2)+64*nf+128*nc+65536;
initialEstimate=nb*272+2*factorBound+n*128+65536;
assert(factorBound<=options.maximumFFFactorBytes&&initialEstimate<=options.maximumConstructionBytes,'SpatialReduction:Budget','原C分块和F因子构造超预算。');
owners=repelem((1:n)',diff(double(packed.rowPointer)));rowCV=packed.order(owners);colCV=packed.order(packed.column);
rowC=split.coarseIndex(rowCV);colC=split.coarseIndex(colCV);rowF=split.fineIndex(rowCV);colF=split.fineIndex(colCV);
B=pick(rowF>0&colC>0,rowF,colC,nf,nc);E=pick(rowC>0&colF>0,rowC,colF,nc,nf);D=pick(rowC>0&colC>0,rowC,colC,nc,nc);
incomingGroup=double(split.fineGroupID(E.column));outgoingGroup=double(split.fineGroupID(B.owner));
CF=spones(sparse(double(E.owner),incomingGroup,ones(numel(E.owner),1),nc,nc));FC=spones(sparse(outgoingGroup,double(B.column),ones(numel(B.column),1),nc,nc));
pairPaths=full(sum(sum(spones(CF),1)'.*sum(spones(FC),2)));
symbolicUpper=pairPaths+numel(D.column);prePatternEstimate=initialEstimate+symbolicUpper*256+nb*48;
assert(prePatternEstimate<=options.maximumConstructionBytes,'SpatialReduction:Budget','所有Schur路径的构造数组超预算。');
pattern=spones(CF*FC+sparse(double(D.owner),double(D.column),ones(numel(D.owner),1),nc,nc));[rr,cc]=find(pattern);pairs=sortrows([rr cc]);ns=size(pairs,1);
assert(ns<=options.maximumSchurBlocks,'SpatialReduction:Budget','完整Schur结构超预算，不丢填充。');
storageBound=factorBound+216*(numel(B.column)+numel(E.column)+numel(D.column))+216*ns+n*96+nc*128+65536;
constructionBound=prePatternEstimate+2*storageBound;
assert(storageBound<=options.maximumReductionBytes&&constructionBound<=options.maximumConstructionBytes,'SpatialReduction:Budget','全部子块/因子/映射或构造数组超预算。');
lookup=sparse(pairs(:,1),pairs(:,2),(1:ns)',nc,nc);values=zeros(5,5,ns);
positions=full(lookup(sub2ind([nc nc],double(D.owner),double(D.column))));
for a=1:5,for b=1:5,values(a,b,positions)=D.values(a,b,:);end,end
fopt=struct('maximumFactorBytes',options.maximumFFFactorBytes,'maximumConstructionBytes',options.maximumConstructionBytes,'failureEvidenceFile',options.failureEvidenceFile);
Ffactor=factorRocketSpatialFBlocks(packed,split,fopt);
[~,inOrder]=sort(incomingGroup);inPointer=[1;1+cumsum(accumarray(incomingGroup,1,[nc 1]))];
[~,outOrder]=sort(outgoingGroup);outPointer=[1;1+cumsum(accumarray(outgoingGroup,1,[nc 1]))];
for g=1:nc
    k=fCounts(g);if k==0,continue;end
    ie=inOrder(inPointer(g):inPointer(g+1)-1);oe=outOrder(outPointer(g):outPointer(g+1)-1);if isempty(ie)||isempty(oe),continue;end
    [coarseRows,~,rowPosition]=unique(double(E.owner(ie)));[coarseCols,~,colPosition]=unique(double(B.column(oe)));
    left=zeros(5*numel(coarseRows),5*k);right=zeros(5*k,5*numel(coarseCols));
    fineCol=double(split.finePosition(E.column(ie)));fineRow=double(split.finePosition(B.owner(oe)));
    for a=1:5,for b=1:5
        idx=sub2ind(size(left),5*(rowPosition-1)+a,5*(fineCol-1)+b);left(idx)=reshape(E.values(a,b,ie),[],1);
        idx=sub2ind(size(right),5*(fineRow-1)+a,5*(colPosition-1)+b);right(idx)=reshape(B.values(a,b,oe),[],1);
    end,end
    inverse=Ffactor.buckets{k}.inverse(:,:,Ffactor.pageOfGroup(g));update=-(left*inverse)*right;
    [columnGrid,rowGrid]=meshgrid(coarseCols,coarseRows);positions=full(lookup(sub2ind([nc nc],rowGrid(:),columnGrid(:))));
    assert(all(positions>0),'SpatialReduction:Pattern','真实局部消元位置缺失，不能删填充。');
    for a=1:5,for b=1:5,values(a,b,positions)=values(a,b,positions)+reshape(update(a:5:end,b:5:end),1,1,[]);end,end
end
assert(all(isfinite(values),'all'),'SpatialReduction:Nonfinite','完整近似Schur系数非有限。');
diagonal=find(pairs(:,1)==pairs(:,2));assert(numel(diagonal)==nc,'SpatialReduction:Diagonal','S缺少原C点对角位置。');
schur=struct('kind','complete-five-conservation-block-CSR-unshifted-steady-graph','cellCount',nc, ...
    'rowPointer',uint32([1;1+cumsum(accumarray(pairs(:,1),1,[nc 1]))]),'column',uint32(pairs(:,2)), ...
    'diagonalPosition',uint32(diagonal),'order',(1:nc)','inverseOrder',(1:nc)','values',values);
reduction=struct('kind','five-field-spatial-local-F-approximate-Schur-reduction','cellCount',n,'split',split,'Ffactor',Ffactor, ...
    'B',B,'E',E,'D',D,'schur',schur,'options',options, ...
    'constructionReport',struct('coarseCells',nc,'FCells',nf,'originalBlocks',nb,'CFBlocks',numel(E.column),'FCBlocks',numel(B.column), ...
    'CCBlocks',numel(D.column),'allLocalEliminationPairPaths',pairPaths,'schurBlocks',ns, ...
    'storageUpperBoundBytes',storageBound,'constructionArrayEstimateBytes',constructionBound,'processPeakMemoryClaimed',false, ...
    'schurIdentity','D-E*M_F*B-with-certified-within-original-patch-F-map', ...
    'fullFFInverseOrExactSchurClaimed',false,'crossPatchFFKeptInCompleteOuterC',true,'smallOrZeroCoefficientsDropped',false, ...
    'allOriginalCVsAndFiveFieldsKept',true,'physicalGridChanged',false,'coarseInverseConstructed',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0));
storage=whos('reduction');assert(storage.bytes<=storageBound&&storage.bytes<=options.maximumReductionBytes,'SpatialReduction:Budget','实际完整空间消元对象超过预算或上界。');
fprintf('SPATIAL_REDUCTION F=%d C=%d SchurBlocks=%d bytes=%.0f\n',nf,nc,ns,storage.bytes);
    function coupling=pick(mask,rows,columns,nrows,ncols)
        coupling=struct('kind','all-original-five-field-rectangular-coupling','rowCount',nrows,'columnCount',ncols, ...
            'owner',uint32(rows(mask)),'column',uint32(columns(mask)),'values',packed.values(:,:,mask));
    end
end
