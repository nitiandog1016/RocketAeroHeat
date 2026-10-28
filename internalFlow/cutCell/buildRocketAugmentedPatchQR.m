function factor=buildRocketAugmentedPatchQR(packed,partitionFactor,metric,penalty)
%BUILDROCKETAUGMENTEDPATCHQR 对完整增广像作局部分块QR，改善相关方向的数值尺度。
% 模块一：沿用原强连接分组，收集每组所有原C行，包含组外行；不删除跨组耦合。
% 模块二：为[C(:,组);penalty*W(:,组)]先逐列归一再作薄QR，以diag(1/norm)*R逆作右映射。
% 模块三：逐组核对增广像正交性。该映射只改变线性求解坐标，不改目标、方程或物理门槛。
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    strcmp(partitionFactor.kind,'fixed-five-conservation-strong-connection-local-block-inverse')&& ...
    packed.cellCount==partitionFactor.cellCount&&metric.cellCount==packed.cellCount, ...
    'AugmentedQR:Input','需要同一当前系统的完整C、原分组及状态度量。');
validateattributes(penalty,{'numeric'},{'scalar','real','finite','positive'});
n=packed.cellCount;nb=numel(packed.column);K=numel(partitionFactor.buckets);
counts=accumarray(double(packed.column),1,[n 1]);
workspaceBound=32*nb+128*n+32*25*K^2*max(counts)+65536;
assert(workspaceBound<=256*1024^2,'AugmentedQR:Budget','局部QR额外构造数组超过256MiB声明上界。');
[~,entries]=sort(packed.column);pointer=[1;1+cumsum(counts)];
owners=repelem((1:n)',diff(double(packed.rowPointer)));buckets=cell(K,1);
coverage=zeros(n,1);used=0;groups=0;maximumDefect=0;minimumQRDiagonal=1;maximumPanelBytes=0;timer=tic;
for k=1:K
    old=partitionFactor.buckets{k};count=numel(old.patchIDs);
    bucket=struct('patchIDs',old.patchIDs,'members',old.members,'map',zeros(5*k,5*k,count), ...
        'gramDefect',zeros(count,1));
    for page=1:count
        physical=double(old.members(:,page));algebraic=double(packed.inverseOrder(physical));
        selected=[];localColumn=[];
        for j=1:k
            ids=double(entries(pointer(algebraic(j)):pointer(algebraic(j)+1)-1));
            selected=[selected;ids(:)];localColumn=[localColumn;repmat(j,numel(ids),1)]; %#ok<AGROW>
        end
        [rows,~,localRow]=unique(owners(selected));panel=zeros(5*numel(rows),5*k);
        for a=1:5,for b=1:5
            index=sub2ind(size(panel),5*(localRow-1)+a,5*(localColumn-1)+b);
            panel(index)=reshape(packed.values(a,b,selected),[],1);
        end,end
        bottom=zeros(5*k);
        for j=1:k
            cellID=physical(j);W=diag(metric.diagonal(cellID,:));
            W(2:4,1)=metric.densityToVelocity(cellID,:)';W(5,1)=metric.densityToTemperature(cellID);
            W(5,2:4)=metric.momentumToTemperature(cellID,:);
            ids=5*(j-1)+(1:5);bottom(ids,ids)=penalty*W;
        end
        augmented=[panel;bottom];norms=vecnorm(augmented,2,1);
        assert(all(isfinite(norms)&norms>0),'AugmentedQR:Column','完整增广像存在零列或非有限列。');
        [~,R]=qr(augmented./norms,0);
        assert(all(isfinite(R),'all')&&all(diag(R)~=0),'AugmentedQR:Rank','局部QR未得到有限非零主元。');
        mapping=(R\eye(5*k))./norms';image=augmented*mapping;
        defect=norm(image'*image-eye(5*k),Inf);
        assert(isfinite(defect)&&defect<=1e-8,'AugmentedQR:Certificate','局部增广像未过原1e-8算术核对门槛。');
        assert(all(coverage(physical)==0),'AugmentedQR:Coverage','原CV被分组重复覆盖。');
        bucket.map(:,:,page)=mapping;bucket.gramDefect(page)=defect;coverage(physical)=1;
        used=used+numel(selected);groups=groups+1;maximumDefect=max(maximumDefect,defect);
        minimumQRDiagonal=min(minimumQRDiagonal,min(abs(diag(R))));
        storage=whos('panel','bottom','augmented','mapping','image');maximumPanelBytes=max(maximumPanelBytes,sum([storage.bytes]));
    end
    buckets{k}=bucket;
end
assert(all(coverage==1)&&used==nb,'AugmentedQR:Coverage','完整原C块或CV未恰好覆盖一次。');
report=struct('method','complete-augmented-image-patch-column-normalization-thin-QR-right-map', ...
    'originalBlocksUsedExactlyOnce',used,'patchCount',groups,'allOriginalCVs',n,'maximumPatchCells',K, ...
    'maximumAugmentedLocalGramDefect',maximumDefect,'minimumNormalizedQRDiagonal',minimumQRDiagonal, ...
    'maximumDensePanelArrayBytes',maximumPanelBytes,'additionalConstructionArrayBoundBytes',workspaceBound, ...
    'normalEquationsFormed',false,'rankTruncationUsed',false,'diagonalCompensationApplied',false, ...
    'globalConditionNumberClaimed',false,'processPeakMemoryClaimed',false,'constructionSeconds',toc(timer), ...
    'originalParentResidualNormUnchanged',true,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
factor=struct('kind','fixed-complete-augmented-image-patch-QR-right-map','cellCount',n,'buckets',{buckets}, ...
    'penaltyInOriginalResidualUnits',penalty,'constructionReport',report);
storage=whos('factor');assert(storage.bytes<=512*1024^2,'AugmentedQR:Budget','持久QR映射超过512MiB上界。');
fprintf('AUGMENTED_QR cells=%d groups=%d blocks=%d maxGram=%.9g\n',n,groups,used,maximumDefect);
end
