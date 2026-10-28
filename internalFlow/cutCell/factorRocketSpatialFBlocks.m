function factor=factorRocketSpatialFBlocks(packed,split,options)
%FACTORROCKETSPATIALFBLOCKS 认证原组内F-F五量主块逆，仅作完整F系统的局部近似。
% 模块一：移除C代表后F组可不连通，不冒称新的强连接物理网格。
% 模块二：局部原逆不合格时显式采用认证的原CV五量对角预条件；阈值1e-8不变。
% 模块三：跨原组F-F仍在完整外层C中，局部因子不能称完整A_FF的逆。
if nargin<3,options=struct;end
defaults=struct('maximumFactorBytes',384*1024^2,'maximumConstructionBytes',3072*1024^2,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'SpatialF:Options','未知F块选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
n=split.cellCount;nf=split.fineCells;ng=split.coarseCells;counts=double(split.groupFCounts);K=max([0;counts]);nb=numel(packed.column);
bound=8*25*sum(counts.^2)+64*nf+128*ng+65536;construction=nb*96+2*bound+64*n;
assert(bound<=options.maximumFactorBytes&&construction<=options.maximumConstructionBytes,'SpatialF:Budget','全部局部F块超预算，不删原CV。');
buckets=cell(K,1);pageOfGroup=zeros(ng,1,'uint32');
for k=1:K
    groups=find(counts==k);pageOfGroup(groups)=uint32(1:numel(groups));members=zeros(k,numel(groups),'uint32');
    for page=1:numel(groups),g=groups(page);members(:,page)=uint32(double(split.groupFPointer(g)):double(split.groupFPointer(g+1))-1);end
    buckets{k}=struct('groupIDs',uint32(groups),'members',members,'inverse',zeros(5*k,5*k,numel(groups)), ...
        'certifiedBounds',zeros(numel(groups),1),'ordinaryDefects',zeros(numel(groups),1), ...
        'fallback',false(numel(groups),1),'rejectedFullFBounds',nan(numel(groups),1));
end
owners=repelem((1:n)',diff(double(packed.rowPointer)));rowCV=packed.order(owners);colCV=packed.order(packed.column);
rowF=split.fineIndex(rowCV);colF=split.fineIndex(colCV);keep=rowF>0&colF>0&split.originalGroupID(rowCV)==split.originalGroupID(colCV);selected=find(keep);
for k=1:K
    ids=selected(counts(split.originalGroupID(rowCV(selected)))==k);if isempty(ids),continue;end
    g=split.originalGroupID(rowCV(ids));page=double(pageOfGroup(g));rr=double(split.finePosition(rowF(ids)));cc=double(split.finePosition(colF(ids)));N=5*k;pages=numel(buckets{k}.groupIDs);
    for a=1:5,for b=1:5,index=sub2ind([N N pages],5*(rr-1)+a,5*(cc-1)+b,page);buckets{k}.inverse(index)=reshape(packed.values(a,b,ids),[],1);end,end
end
included=numel(selected);clear owners rowCV colCV rowF colF keep selected ids page rr cc index;
maximum=0;fallbacks=0;accurate=0;refined=0;
for k=1:K
    bucket=buckets{k};buckets{k}=struct([]);
    for page=1:numel(bucket.groupIDs)
        original=bucket.inverse(:,:,page);[candidate,certificate]=certifyRocketLocalInverse(original);
        accurate=accurate+certificate.accurateResidualUsed;refined=refined+(certificate.refinementSteps>0);
        if ~certificate.usable
            rejected=certificate.upperBound;replacement=zeros(5*k);upper=0;ordinary=0;
            for member=1:k
                ids=5*(member-1)+(1:5);[local,localCertificate]=certifyRocketLocalInverse(original(ids,ids));
                if ~localCertificate.usable
                    groupID=bucket.groupIDs(page);originalCV=split.fineOriginalCVs(bucket.members(member,page));
                    if ~isempty(options.failureEvidenceFile),save(options.failureEvidenceFile,'original','groupID','originalCV','candidate','certificate','localCertificate','-v7');end
                    error('SpatialF:UnsafeCV','原F控制体%d五量对角逆未通过1e-8。',originalCV);
                end
                replacement(ids,ids)=local;upper=max(upper,localCertificate.upperBound);ordinary=max(ordinary,localCertificate.ordinaryBLASDefect);
            end
            candidate=replacement;certificate.upperBound=upper;certificate.ordinaryBLASDefect=ordinary;
            bucket.fallback(page)=true;bucket.rejectedFullFBounds(page)=rejected;fallbacks=fallbacks+1;
        end
        bucket.inverse(:,:,page)=candidate;bucket.certifiedBounds(page)=certificate.upperBound;bucket.ordinaryDefects(page)=certificate.ordinaryBLASDefect;maximum=max(maximum,certificate.upperBound);
    end
    buckets{k}=bucket;clear bucket;
end
factor=struct('kind','fixed-within-original-patch-five-field-F-inverse','cellCount',nf,'groupCount',ng,'buckets',{buckets},'pageOfGroup',pageOfGroup,'groupFCounts',split.groupFCounts, ...
    'factorizationReport',struct('includedOriginalFFBlocks',included,'maximumCertifiedOriginalInputInverseDefect',maximum,'inverseDefectThreshold',1e-8, ...
    'fallbackFBlockCount',fallbacks,'accurateCheckedFBlocks',accurate,'refinedFBlocks',refined, ...
    'factorStorageUpperBoundBytes',bound,'constructionArrayEstimateBytes',construction,'processPeakMemoryClaimed',false, ...
    'fixedLinearMap',true,'fullFFInverseClaimed',false,'rejectedFullFInverseClaimedValid',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0));
storage=whos('factor');assert(storage.bytes<=bound&&storage.bytes<=options.maximumFactorBytes,'SpatialF:Budget','实际F因子超过预算或上界。');
fprintf('SPATIAL_F cells=%d groups=%d fallback=%d bound=%.9g bytes=%.0f\n',nf,ng,fallbacks,maximum,storage.bytes);
end
