function report=testRocketFullFStrength(folder)
%TESTROCKETFULLFSTRENGTH 两排序核对F图重新分组的固定局部逆及完整F/原C参考。
% 模块一：只在完整A_FF上建立代数组，原F/C角色及物理网格不变。
% 模块二：全部五量缓存逆、完整F方向和固定C回代分别检查，不沿用旧S身份。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FullFStrengthTest:Overwrite','不得覆盖F强连接小域。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','fullFunit01','full_F_unit_evidence.mat');s=load(parentFile,'cases');cases=cell(2,1);
for j=1:2
    source=s.cases{j};graph=source.graph;before=graph;partition=buildRocketStrongPatches(graph,4);
    factor=factorRocketCoupledPatches(graph,partition,struct('maximumFactorBytes',8*1024^2,'maximumConstructionBytes',64*1024^2));
    assert(isequaln(before,graph));nf=graph.cellCount;P=zeros(5*nf);
    coverage=rocketFullFStrengthCoverage(graph,factor,source.split.fineGroupID);
    coverageRow=rocketFullFStrengthCoverage(graph,factor,source.split.fineGroupID');
    assert(isequaln(coverage,coverageRow)&&coverage.oldCrossPatchBlocks==source.extraction.crossPatchFFBlocks);
    for k=1:numel(factor.buckets)
        bucket=factor.buckets{k};for page=1:numel(bucket.patchIDs)
            members=double(bucket.members(:,page));ids=reshape(members(:)'+(0:4)'*nf,[],1);block=source.A(ids,ids);
            if bucket.fallback(page),diagonal=zeros(5*k);for r=1:k,ii=5*(r-1)+(1:5);diagonal(ii,ii)=block(ii,ii);end;block=diagonal;end
            inverse=block\eye(5*k);assert(norm(inverse-bucket.inverse(:,:,page),'fro')/max(1,norm(inverse,'fro'))<1e-8);P(ids,ids)=inverse;
        end
    end
    rhs=source.rhs;single=solveRocketCoupledPatches(factor,rhs);assert(norm(single(:)-P*rhs(:))/max(1,norm(P*rhs(:)))<1e-9);
    opt=struct('maximumIterations',5*nf,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
    [value,detail]=solveRocketCoupledPatchKrylov(graph,factor,rhs,opt);assert(detail.krylov.flag==0&&detail.relativeCResidual<1.1e-10);
    direct=source.A\rhs(:);assert(norm(value(:)-direct)/max(1,norm(direct))<1e-8);
    complete=source.complete;complete(source.split.fineOriginalCVs,:)=value;parentImage=applyRocketPackedGraph(source.packed,complete);
    parentResidual=norm(source.parentRhs-parentImage,'fro')/norm(source.parentRhs,'fro');assert(parentResidual<1e-8);
    cases{j}=struct('graph',graph,'factor',factor,'coverage',coverage,'P',P,'A',source.A,'rhs',rhs,'single',single,'value',value,'detail',detail, ...
        'packed',source.packed,'split',source.split,'complete',complete,'parentRhs',source.parentRhs,'parentImage',parentImage,'parentResidual',parentResidual);
end
report=struct('status','passed-two-orderings-F-only-strength-local-inverse-and-complete-F-parent-reference', ...
    'rowAndColumnOldGroupCoverageExactlyEqual',true,'physicalGridChanged',false,'oldApproximateSchurIdentityClaimed',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'full_F_strength_unit_evidence.mat'),'report','cases','parentFile','-v7');fid=fopen(fullfile(folder,'full_F_strength_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
