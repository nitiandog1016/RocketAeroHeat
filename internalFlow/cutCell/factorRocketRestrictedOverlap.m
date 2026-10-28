function factor=factorRocketRestrictedOverlap(graph,input,options)
%FACTORROCKETRESTRICTEDOVERLAP 认证固定扩展F主块逆，仅将原核心行用于预条件输出。
% 模块一：核心点恰一次、核心成员在扩展组前缀，至多两个邻点，原物理网格不变。
% 模块二：按核心/扩展大小分桶，每组从完整F图取全部原系数，不丢小非零。
% 模块三：原扩展逆未过1e-8时显式改用认证原CV五量对角，旧整组逆仍拒绝。
if nargin<3,options=struct;end
defaults=struct('maximumFactorBytes',512*1024^2,'maximumConstructionBytes',1536*1024^2,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'OverlapF:Options','未知受限重叠选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
for key={'maximumFactorBytes','maximumConstructionBytes'},validateattributes(options.(key{1}),{'numeric'},{'scalar','real','finite','positive'});end
assert(strcmp(graph.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    isfield(graph,'systemScope')&&strcmp(graph.systemScope,'exact-spatial-F-submatrix-A_FF-not-original-full-C'),'OverlapF:Input','需要已验完整A_FF图。');
required={'coreGroupIDOneBased','coreMembersOneBased','corePointerOneBased','overlapMembersOneBased','overlapPointerOneBased','coreSizes','overlapSizes'};
assert(isstruct(input)&&isscalar(input)&&all(isfield(input,required)),'OverlapF:Input','核心/扩展成员输入不完整。');
t=struct;
for key=required
    raw=input.(key{1});assert(isnumeric(raw)&&isreal(raw)&&isvector(raw)&&all(isfinite(raw))&&all(raw>=1&raw==floor(raw)),'OverlapF:Input','成员或指针必须为有限正整数向量。');
    t.(key{1})=double(raw(:));
end
n=graph.cellCount;ng=numel(t.coreSizes);m=t.coreSizes;k=t.overlapSizes;K=max(k);M=max(m);nb=numel(graph.column);
assert(numel(t.coreGroupIDOneBased)==n&&numel(t.coreMembersOneBased)==n&&isequal(sort(t.coreMembersOneBased),(1:n)')&& ...
    numel(t.corePointerOneBased)==ng+1&&t.corePointerOneBased(1)==1&&t.corePointerOneBased(end)==n+1&& ...
    isequal(diff(t.corePointerOneBased),m)&&all(m<=8)&&all(t.coreGroupIDOneBased<=ng),'OverlapF:Coverage','原核心控制体必须恰好覆盖所有F点一次。');
assert(numel(k)==ng&&all(k>=m&k<=m+2)&&numel(t.overlapPointerOneBased)==ng+1&& ...
    t.overlapPointerOneBased(1)==1&&t.overlapPointerOneBased(end)==numel(t.overlapMembersOneBased)+1&& ...
    isequal(diff(t.overlapPointerOneBased),k)&&all(t.overlapMembersOneBased<=n),'OverlapF:Members','扩展组大小或成员范围非法。');
for g=1:ng
    core=t.coreMembersOneBased(t.corePointerOneBased(g):t.corePointerOneBased(g+1)-1);
    members=t.overlapMembersOneBased(t.overlapPointerOneBased(g):t.overlapPointerOneBased(g+1)-1);
    assert(all(t.coreGroupIDOneBased(core)==g)&&numel(unique(members))==k(g)&&isequal(members(1:m(g)),core),'OverlapF:Members','每组核心必须为扩展成员的唯一前缀。');
end
memberCount=numel(t.overlapMembersOneBased);factorBound=8*25*sum(k.^2)+96*n+160*ng+8*memberCount+65536;
construction=200*nb+2*factorBound+128*n+48*memberCount+256*ng+65536;
assert(factorBound<=options.maximumFactorBytes&&construction<=options.maximumConstructionBytes,'OverlapF:Budget','全部扩展逆或构造数组超预算，不能删核心/邻点。');
topology=struct;for key=required
    if ismember(key{1},{'coreSizes','overlapSizes'}),topology.(key{1})=uint8(t.(key{1}));else,topology.(key{1})=uint32(t.(key{1}));end
end
bucketID=(m-1)*K+k;buckets=cell(M*K,1);
for id=1:M*K
    groups=find(bucketID==id);coreSize=floor((id-1)/K)+1;overlapSize=mod(id-1,K)+1;
    members=zeros(overlapSize,numel(groups),'uint32');cores=zeros(coreSize,numel(groups),'uint32');
    for page=1:numel(groups)
        g=groups(page);members(:,page)=topology.overlapMembersOneBased(double(topology.overlapPointerOneBased(g)):double(topology.overlapPointerOneBased(g+1))-1);
        cores(:,page)=topology.coreMembersOneBased(double(topology.corePointerOneBased(g)):double(topology.corePointerOneBased(g+1))-1);
    end
    buckets{id}=struct('coreSize',coreSize,'overlapSize',overlapSize,'groupIDs',uint32(groups),'members',members,'coreMembers',cores, ...
        'inverse',zeros(5*overlapSize,5*overlapSize,numel(groups)),'certifiedBounds',zeros(numel(groups),1), ...
        'ordinaryDefects',zeros(numel(groups),1),'fallback',false(numel(groups),1),'rejectedFullBounds',nan(numel(groups),1));
end
included=0;maximum=0;accurate=0;refined=0;fallbacks=0;cvChecks=0;timer=tic;
for id=1:numel(buckets)
    bucket=buckets{id};buckets{id}=struct([]);extent=bucket.overlapSize;
    for page=1:numel(bucket.groupIDs)
        members=double(bucket.members(:,page));original=zeros(5*extent);
        for localRow=1:extent
            row=graph.inverseOrder(members(localRow));edges=double(graph.rowPointer(row)):double(graph.rowPointer(row+1))-1;
            physicalColumns=graph.order(graph.column(edges));[keep,localColumns]=ismember(physicalColumns,members);edges=edges(keep);localColumns=localColumns(keep);
            included=included+numel(edges);
            for e=1:numel(edges),rr=5*(localRow-1)+(1:5);cc=5*(localColumns(e)-1)+(1:5);original(rr,cc)=graph.values(:,:,edges(e));end
        end
        [candidate,certificate]=certifyRocketLocalInverse(original);accurate=accurate+certificate.accurateResidualUsed;refined=refined+(certificate.refinementSteps>0);
        if ~certificate.usable
            rejectedBound=certificate.upperBound;replacement=zeros(5*extent);upper=0;ordinary=0;
            for member=1:extent
                ii=5*(member-1)+(1:5);[local,localCertificate]=certifyRocketLocalInverse(original(ii,ii));cvChecks=cvChecks+1;
                if ~localCertificate.usable
                    groupID=bucket.groupIDs(page);originalFPoint=members(member);
                    if ~isempty(options.failureEvidenceFile),save(options.failureEvidenceFile,'original','members','groupID','originalFPoint','candidate','certificate','localCertificate','-v7');end
                    error('OverlapF:UnsafeCV','F点%d的五量对角逆未过1e-8，整重叠映射拒绝。',originalFPoint);
                end
                replacement(ii,ii)=local;upper=max(upper,localCertificate.upperBound);ordinary=max(ordinary,localCertificate.ordinaryBLASDefect);
            end
            candidate=replacement;certificate.upperBound=upper;certificate.ordinaryBLASDefect=ordinary;
            bucket.fallback(page)=true;bucket.rejectedFullBounds(page)=rejectedBound;fallbacks=fallbacks+1;
        end
        bucket.inverse(:,:,page)=candidate;bucket.certifiedBounds(page)=certificate.upperBound;bucket.ordinaryDefects(page)=certificate.ordinaryBLASDefect;
        maximum=max(maximum,certificate.upperBound);
    end
    buckets{id}=bucket;clear bucket;
end
factor=struct('kind','fixed-five-field-overlap-inverse-return-original-core-rows-only','cellCount',n,'groupCount',ng,'topology',topology,'buckets',{buckets},'options',options, ...
    'factorizationReport',struct('corePatches',ng,'maximumCorePoints',M,'maximumOverlapPoints',K,'overlapMembers',memberCount, ...
    'originalFFBlockOccurrencesInsideOverlapSets',included,'inverseNumericArrayBytes',8*25*sum(k.^2), ...
    'factorStorageUpperBoundBytes',factorBound,'constructionArrayEstimateBytes',construction,'processPeakMemoryClaimed',false, ...
    'maximumCertifiedOriginalInputInverseDefect',maximum,'inverseDefectThreshold',1e-8,'fallbackOverlapBlocks',fallbacks, ...
    'accurateCheckedFullBlocks',accurate,'refinedFullBlocks',refined,'fallbackCVDiagonalInverseChecks',cvChecks,'factorSeconds',toc(timer), ...
    'returnOnlyOriginalCoreRows',true,'coreUnknownsReturnedExactlyOnce',true,'fixedLinearMap',true,'physicalGridChanged',false, ...
    'rejectedFullOverlapInverseClaimedValid',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0));
storage=whos('factor');assert(storage.bytes<=factorBound&&storage.bytes<=options.maximumFactorBytes,'OverlapF:Budget','实际重叠缓存超过预算或上界。');
fprintf('RESTRICTED_OVERLAP cells=%d groups=%d max=%d fallback=%d bound=%.9g bytes=%.0f\n',n,ng,K,fallbacks,maximum,storage.bytes);
end
