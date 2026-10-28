function hierarchy=buildRocketPatchHierarchy(packed,rootFactor,options)
%BUILDROCKETPATCHHIERARCHY 复用已验五量局部逆，建立保留全部系数的强连接代数层次。
% 模块一：每层完整C、固定局部逆及R/P均保留；只在预条件中使用粗系统。
% 模块二：最粗完整五量系统必须通过原输入逆残差1e-8，失败不补对角。
% 模块三：同时限制层数、全部层存储及构造数组；不缩物理域或放宽原方程。
if nargin<3,options=struct;end
defaults=struct('maximumPatchCells',8,'maximumLevels',8,'terminalMaximumCells',16, ...
    'maximumHierarchyBytes',1024*1024^2,'maximumConstructionBytes',3072*1024^2,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'PatchHierarchy:Options','未知代数层次选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumPatchCells,{'numeric'},{'scalar','integer','>=',2,'<=',16});
validateattributes(options.maximumLevels,{'numeric'},{'scalar','integer','>=',2,'<=',12});
validateattributes(options.terminalMaximumCells,{'numeric'},{'scalar','integer','>=',1,'<=',32});
for key={'maximumHierarchyBytes','maximumConstructionBytes'},validateattributes(options.(key{1}),{'numeric'},{'scalar','finite','real','positive'});end
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    strcmp(rootFactor.kind,'fixed-five-conservation-strong-connection-local-block-inverse')&&rootFactor.cellCount==packed.cellCount,'PatchHierarchy:Input','需要同一完整C及已验固定五量局部逆。');
levels=cell(options.maximumLevels,1);count=0;current=packed;rootCells=packed.cellCount;newFactors=0;transfers=struct([]);
while true
    count=count+1;assert(count<=options.maximumLevels,'PatchHierarchy:Levels','层数预算用尽，不能假装最粗已解。');
    n=current.cellCount;
    if n<=options.terminalMaximumCells
        N=5*n;dense=zeros(N);
        for row=1:n
            physicalRow=current.order(row);ids=double(current.rowPointer(row)):double(current.rowPointer(row+1))-1;
            for edge=ids,physicalColumn=current.order(current.column(edge));dense(physicalRow+(0:4)*n,physicalColumn+(0:4)*n)=current.values(:,:,edge);end
        end
        [inverse,certificate]=certifyRocketLocalInverse(dense);
        if ~certificate.usable
            if ~isempty(options.failureEvidenceFile),save(options.failureEvidenceFile,'dense','inverse','certificate','count','n','-v7');end
            error('PatchHierarchy:TerminalInverse','最粗完整五量逆未通过原输入1e-8，拒绝整层次。');
        end
        levels{count}=struct('terminal',true,'packed',current,'factor',struct([]),'inverse',inverse,'certificate',certificate);
        break;
    end
    if count==1
        factor=rootFactor;
    else
        partition=buildRocketStrongPatches(current,options.maximumPatchCells);
        factorOptions=struct('maximumFactorBytes',options.maximumHierarchyBytes,'maximumConstructionBytes',options.maximumConstructionBytes, ...
            'failureEvidenceFile',options.failureEvidenceFile);
        factor=factorRocketCoupledPatches(current,partition,factorOptions);newFactors=newFactors+1;
    end
    assert(factor.partition.patchCount<n,'PatchHierarchy:Coarsening','没有减少代数未知量，拒绝循环。');
    levels{count}=struct('terminal',false,'packed',current,'factor',factor,'inverse',[],'certificate',struct([]));
    live=whos('levels');assert(live.bytes<=options.maximumHierarchyBytes,'PatchHierarchy:Budget','完整层存储超预算。');
    [coarse,report]=coarsenRocketPackedGraph(current,factor.partition,options.maximumConstructionBytes);
    probe=reshape(sin((1:5*coarse.cellCount)'*.317),coarse.cellCount,5);
    lifted=prolongRocketPatchCorrection(factor.partition,probe);roundtrip=restrictRocketPatchResidual(factor.partition,lifted);
    fineImage=applyRocketPackedGraph(current,lifted);restricted=restrictRocketPatchResidual(factor.partition,fineImage);coarseImage=applyRocketPackedGraph(coarse,probe);
    report.RPRelativeDifference=norm(roundtrip-probe,'fro')/max(1,norm(probe,'fro'));
    report.RCPActionRelativeDifference=norm(restricted-coarseImage,'fro')/max(1,norm(coarseImage,'fro'));
    assert(report.RPRelativeDifference<1e-12&&report.RCPActionRelativeDifference<1e-8,'PatchHierarchy:Transfer','完整R/P或原系数聚合作用不符。');
    report.fineCVerificationApplications=1;report.coarseCVerificationApplications=1;
    if isempty(transfers),transfers=report;else,transfers(end+1)=report;end %#ok<AGROW>
    current=coarse;
end
levels=levels(1:count);cellCounts=zeros(count,1);blockCounts=cellCounts;
for j=1:count,cellCounts(j)=levels{j}.packed.cellCount;blockCounts(j)=numel(levels{j}.packed.column);end
hierarchy=struct('kind','fixed-all-five-field-unsmoothed-strong-connection-V-cycle','cellCount',rootCells,'levels',{levels},'options',options, ...
    'constructionReport',struct('levels',count,'cellCounts',cellCounts,'blockCounts',blockCounts,'newLocalFactors',newFactors, ...
    'rootFactorReused',true,'transfers',transfers,'terminalInverseCertificate',certificate, ...
    'cycle','one-fixed-local-pre-smooth-recursive-coarse-correction-no-post-smooth', ...
    'fixedLinearMap',true,'allOriginalPhysicalCVsAndFiveFieldsKept',true,'originalFineCChanged',false, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'processPeakMemoryClaimed',false));
storage=whos('hierarchy');assert(storage.bytes<=options.maximumHierarchyBytes,'PatchHierarchy:Budget','实际完整层次超过存储预算。');
fprintf('PATCH_HIERARCHY levels=%d cells=%s bytes=%.0f terminalBound=%.9g\n',count,mat2str(cellCounts'),storage.bytes,certificate.upperBound);
end
