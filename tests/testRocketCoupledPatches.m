function report=testRocketCoupledPatches(folder)
%TESTROCKETCOUPLEDPATCHES 两种排列核对强连接组、全部五量缓存逆及完整C稠密参考。
% 模块一：所有原CV恰一次、每组沿原非零连接连通，原C不因局部分组而删项。
% 模块二：逐块直接逆、固定线性映射、完整C和物理Arnoldi方向共同核对。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CoupledPatchesTest:Overwrite','不得覆盖五量局部块测试。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','blockiluunit01','block_ilu_unit_evidence.mat');s=load(parentFile);cases=cell(2,1);
for j=1:2
    packed=s.cases{j}.packed;before=packed;partition=buildRocketStrongPatches(packed,4);n=packed.cellCount;assert(isequal(sort(double(partition.members)),(1:n)'));
    options=struct('maximumFactorBytes',8*1024^2,'maximumConstructionBytes',64*1024^2);factor=factorRocketCoupledPatches(packed,partition,options);assert(isequaln(before,packed));P=zeros(5*n);
    for k=1:numel(factor.buckets)
        bucket=factor.buckets{k};for page=1:numel(bucket.patchIDs)
            member=double(bucket.members(:,page));ids=reshape(member(:)'+(0:4)'*n,[],1);block=s.C(ids,ids);
            if bucket.fallback(page),diagBlock=zeros(5*k);for r=1:k,ii=5*(r-1)+(1:5);diagBlock(ii,ii)=block(ii,ii);end;block=diagBlock;end
            reference=block\eye(5*k);assert(norm(reference-bucket.inverse(:,:,page),'fro')/max(1,norm(reference,'fro'))<1e-8);P(ids,ids)=reference;
        end
    end
    rhs=s.cases{j}.rhs;single=solveRocketCoupledPatches(factor,rhs);assert(norm(single(:)-P*rhs(:))/max(1,norm(P*rhs(:)))<1e-9);
    other=reshape(cos((1:5*n)'*.29),n,5);combined=solveRocketCoupledPatches(factor,.7*rhs-1.3*other);linearity=norm(combined-.7*single+1.3*solveRocketCoupledPatches(factor,other),'fro')/max(1,norm(combined,'fro'));assert(linearity<1e-9);
    opt=struct('maximumIterations',5*n,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);[value,detail]=solveRocketCoupledPatchKrylov(packed,factor,rhs,opt);direct=s.C\rhs(:);
    assert(detail.krylov.flag==0&&detail.relativeCResidual<1.1e-10&&norm(value(:)-direct)/max(1,norm(direct))<1e-8);
    compact=compactRocketCoupledKrylovReport(detail);raw=jsonencode(compact);assert(compact.relativeCResidual==detail.relativeCResidual&&~contains(raw,'"V":')&&~contains(raw,'"Z":')&&numel(raw)<128*1024);
    bad=options;bad.maximumFactorBytes=1;rejected=false;try,factorRocketCoupledPatches(packed,partition,bad);catch ex,rejected=strcmp(ex.identifier,'CoupledPatches:Budget');end;assert(rejected);
    [zero,z]=solveRocketCoupledPatchKrylov(packed,factor,zeros(n,5),opt);assert(all(zero==0,'all')&&z.CApplications==0&&z.localBlockApplications==0);
    cases{j}=struct('packed',packed,'partition',partition,'factor',factor,'P',P,'rhs',rhs,'single',single,'value',value,'detail',detail,'direct',direct,'linearity',linearity);
end
% 模块三：奇异组只在预条件中改用原CV五量对角，完整原C仍含全部耦合。
toy=struct('kind','complete-five-conservation-block-CSR-unshifted-steady-graph','cellCount',2,'rowPointer',uint32([1;3;5]), ...
    'column',uint32([1;2;1;2]),'diagonalPosition',uint32([1;4]),'order',[1;2],'inverseOrder',[1;2],'values',repmat(eye(5),1,1,4));
toyPartition=buildRocketStrongPatches(toy,2);toyFactor=factorRocketCoupledPatches(toy,toyPartition,options);toyRhs=reshape(1:10,2,5);toyValue=solveRocketCoupledPatches(toyFactor,toyRhs);toyImage=applyRocketPackedGraph(toy,toyValue);
assert(toyFactor.factorizationReport.fallbackPatchCount==1&&isequal(toyValue,toyRhs)&&isequal(toyImage,repmat(sum(toyRhs,1),2,1)));
report=struct('status','passed-small-two-orderings-bounded-strong-five-field-patches-and-full-C-reference','explicitCVFallbackKeepsOriginalC',true,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'coupled_patches_unit_evidence.mat'),'report','cases','parentFile','toy','toyPartition','toyFactor','toyRhs','toyValue','toyImage','-v7');fid=fopen(fullfile(folder,'coupled_patches_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
