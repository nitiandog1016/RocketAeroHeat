function report=testRocketRestrictedOverlap(folder)
%TESTROCKETRESTRICTEDOVERLAP 两排序核对扩展块逆、仅核心回写、线性及完整F/原C参考。
% 模块一：小域独立构造重叠组，不使用真实预检来替代数值参考。
% 模块二：零重叠须退化为原强连接局部逆；所有核心控制体恰一次输出。
% 模块三：行/列输入、预算、坏成员、奇异扩展组替代及零对角拒绝单列。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'OverlapFTest:Overwrite','不得覆盖受限重叠测试。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','fullFstrengthunit02','full_F_strength_unit_evidence.mat');s=load(parentFile,'cases');cases=cell(2,1);
options=struct('maximumFactorBytes',32*1024^2,'maximumConstructionBytes',128*1024^2);
for j=1:2
    source=s.cases{j};graph=source.graph;before=graph;topology=makeTopology(graph,source.factor.partition);
    factor=factorRocketRestrictedOverlap(graph,topology,options);assert(isequaln(graph,before));nf=graph.cellCount;P=zeros(5*nf);
    for id=1:numel(factor.buckets)
        bucket=factor.buckets{id};m=bucket.coreSize;k=bucket.overlapSize;
        for page=1:numel(bucket.groupIDs)
            members=double(bucket.members(:,page));core=double(bucket.coreMembers(:,page));ids=reshape(members(:)'+(0:4)'*nf,[],1);rr=reshape(core(:)'+(0:4)'*nf,[],1);
            block=source.A(ids,ids);if bucket.fallback(page),diagonal=zeros(5*k);for c=1:k,ii=5*(c-1)+(1:5);diagonal(ii,ii)=block(ii,ii);end;block=diagonal;end
            inverse=block\eye(5*k);assert(norm(inverse-bucket.inverse(:,:,page),'fro')/max(1,norm(inverse,'fro'))<1e-8);P(rr,ids)=inverse(1:5*m,:);
        end
    end
    rhs=source.rhs;single=solveRocketRestrictedOverlap(factor,rhs);assert(norm(single(:)-P*rhs(:))/max(1,norm(P*rhs(:)))<1e-9);
    other=reshape(cos((1:5*nf)'*.29),nf,5);otherValue=solveRocketRestrictedOverlap(factor,other);combined=solveRocketRestrictedOverlap(factor,.7*rhs-1.3*other);
    assert(norm(combined-.7*single+1.3*otherValue,'fro')/max(1,norm(combined,'fro'))<1e-9);
    rowTopology=topology;for key=fieldnames(rowTopology)',rowTopology.(key{1})=rowTopology.(key{1})';end
    rowFactor=factorRocketRestrictedOverlap(graph,rowTopology,options);assert(isequaln(factor.topology,rowFactor.topology)&&isequal(solveRocketRestrictedOverlap(rowFactor,rhs),single));
    plain=topology;plain.overlapMembersOneBased=plain.coreMembersOneBased;plain.overlapPointerOneBased=plain.corePointerOneBased;plain.overlapSizes=plain.coreSizes;
    plainFactor=factorRocketRestrictedOverlap(graph,plain,options);plainValue=solveRocketRestrictedOverlap(plainFactor,rhs);reference=solveRocketCoupledPatches(source.factor,rhs);
    assert(norm(plainValue-reference,'fro')/max(1,norm(reference,'fro'))<1e-10);
    opt=struct('maximumIterations',5*nf,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
    [value,detail]=solveRocketOverlapFullF(graph,factor,rhs,opt);direct=source.A\rhs(:);
    assert(detail.krylov.flag==0&&detail.relativeFullFResidual<1.1e-10&&norm(value(:)-direct)/max(1,norm(direct))<1e-8);
    complete=source.complete;complete(source.split.fineOriginalCVs,:)=value;parentImage=applyRocketPackedGraph(source.packed,complete);
    parentResidual=norm(source.parentRhs-parentImage,'fro')/norm(source.parentRhs,'fro');assert(parentResidual<1e-8);
    compact=compactRocketFullFReport(detail);raw=jsonencode(compact);assert(~contains(raw,'"V":')&&~contains(raw,'"Z":')&&numel(raw)<128*1024);
    [zero,z]=solveRocketOverlapFullF(graph,factor,zeros(nf,5),opt);assert(all(zero==0,'all')&&z.FApplications==0&&z.localFApplications==0);
    bad=options;bad.maximumFactorBytes=1;rejected=false;try,factorRocketRestrictedOverlap(graph,topology,bad);catch ex,rejected=strcmp(ex.identifier,'OverlapF:Budget');end;assert(rejected);
    bad=topology;bad.coreMembersOneBased(1)=bad.coreMembersOneBased(2);rejected=false;try,factorRocketRestrictedOverlap(graph,bad,options);catch ex,rejected=strcmp(ex.identifier,'OverlapF:Coverage');end;assert(rejected);
    cases{j}=struct('graph',graph,'topology',topology,'factor',factor,'P',P,'A',source.A,'rhs',rhs,'single',single,'value',value,'detail',detail, ...
        'plainFactor',plainFactor,'plainValue',plainValue,'packed',source.packed,'split',source.split,'complete',complete,'parentRhs',source.parentRhs,'parentImage',parentImage,'parentResidual',parentResidual);
end
% 模块四：奇异扩展组只在预条件内退到CV对角，原完整F算子仍保留耦合。
toy=struct('kind','complete-five-conservation-block-CSR-unshifted-steady-graph','cellCount',2,'rowPointer',uint32([1;3;5]), ...
    'column',uint32([1;2;1;2]),'diagonalPosition',uint32([1;4]),'order',[1;2],'inverseOrder',[1;2], ...
    'values',repmat(eye(5),1,1,4),'systemScope','exact-spatial-F-submatrix-A_FF-not-original-full-C');
t=struct('coreGroupIDOneBased',[1;2],'coreMembersOneBased',[1;2],'corePointerOneBased',[1;2;3], ...
    'overlapMembersOneBased',[1;2;2;1],'overlapPointerOneBased',[1;3;5],'coreSizes',[1;1],'overlapSizes',[2;2]);
fallback=factorRocketRestrictedOverlap(toy,t,options);toyRhs=reshape(1:10,2,5);toyValue=solveRocketRestrictedOverlap(fallback,toyRhs);toyImage=applyRocketPackedGraph(toy,toyValue);
assert(fallback.factorizationReport.fallbackOverlapBlocks==2&&isequal(toyValue,toyRhs)&&isequal(toyImage,repmat(sum(toyRhs,1),2,1)));
badToy=toy;badToy.values(:,:,1:2)=0;unsafeRejected=false;try,factorRocketRestrictedOverlap(badToy,t,options);catch ex,unsafeRejected=strcmp(ex.identifier,'OverlapF:UnsafeCV');end;assert(unsafeRejected);
report=struct('status','passed-two-orderings-restricted-overlap-all-inverses-core-output-linearity-and-full-F-reference', ...
    'rowColumnTopologyExactlyEqual',true,'zeroOverlapMatchesPreviousLocalMap',true,'unsafeCVDiagonalRejected',true, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'restricted_overlap_unit_evidence.mat'),'report','cases','parentFile','toy','fallback','toyRhs','toyValue','toyImage','-v7');fid=fopen(fullfile(folder,'restricted_overlap_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function t=makeTopology(graph,partition)
% 模块五：测试用小域扩展按原点编号取最多两个相邻非核心点，所有数组显式列化。
n=graph.cellCount;ng=partition.patchCount;members=[];sizes=zeros(ng,1);pointer=ones(ng+1,1);
for g=1:ng
    core=double(partition.members(double(partition.groupPointer(g)):double(partition.groupPointer(g+1))-1));pool=[];
    for node=core'
        row=graph.inverseOrder(node);ids=double(graph.rowPointer(row)):double(graph.rowPointer(row+1))-1;pool=[pool;double(graph.order(graph.column(ids)))]; %#ok<AGROW>
    end
    pool=setdiff(unique(pool),core);extra=pool(1:min(2,numel(pool)));local=[core(:);extra(:)];members=[members;local];sizes(g)=numel(local);pointer(g+1)=numel(members)+1; %#ok<AGROW>
end
assert(numel(partition.members)==n);
t=struct('coreGroupIDOneBased',double(partition.groupID(:)),'coreMembersOneBased',double(partition.members(:)), ...
    'corePointerOneBased',double(partition.groupPointer(:)),'overlapMembersOneBased',members(:),'overlapPointerOneBased',pointer(:), ...
    'coreSizes',double(partition.sizes(:)),'overlapSizes',sizes(:));
end
