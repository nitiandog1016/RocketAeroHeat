function report=testRocketSpatialReduction(folder)
%TESTROCKETSPATIALREDUCTION 两排序稠密核对空间五量分块、近似Schur和完整回代。
% 模块一：保留所有原C/F控制体，局部F逆与完整A_FF逆明确区分。
% 模块二：全部原矩形块、Schur系数、固定映射和原C的Arnoldi方向共同验收。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SpatialReductionTest:Overwrite','不得覆盖空间消元测试。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','coupledpatchunit03','coupled_patches_unit_evidence.mat');s=load(parentFile,'cases');cases=cell(2,1);
options=struct('maximumReductionBytes',32*1024^2,'maximumFFFactorBytes',16*1024^2,'maximumConstructionBytes',128*1024^2);
for j=1:2
    source=s.cases{j};packed=source.packed;before=packed;reduction=buildRocketSpatialReduction(packed,source.partition,[],options);split=reduction.split;n=split.cellCount;nf=split.fineCells;nc=split.coarseCells;
    assert(isequaln(before,packed));C=densePacked(packed);fids=reshape(double(split.fineOriginalCVs)+(0:4)*n,[],1);cids=reshape(double(split.coarseOriginalCVs)+(0:4)*n,[],1);
    A=C(fids,fids);B=C(fids,cids);E=C(cids,fids);D=C(cids,cids);MF=denseFInverse(reduction.Ffactor,A);
    assert(norm(B-denseRect(reduction.B),'fro')<1e-12&&norm(E-denseRect(reduction.E),'fro')<1e-12&&norm(D-denseRect(reduction.D),'fro')<1e-12);
    S=D-E*MF*B;assembled=densePacked(reduction.schur);assert(norm(assembled-S,'fro')/max(1,norm(S,'fro'))<1e-12);
    probe=reshape(sin((1:5*nc)'*.371),nc,5);implicit=applyRocketSpatialSchur(reduction,probe);assert(norm(implicit(:)-S*probe(:))/max(1,norm(S*probe(:)))<1e-11);
    coarsePartition=buildRocketStrongPatches(reduction.schur,4);coarseFactor=factorRocketCoupledPatches(reduction.schur,coarsePartition, ...
        struct('maximumFactorBytes',8*1024^2,'maximumConstructionBytes',64*1024^2));MS=denseCoarseInverse(coarseFactor,S);
    M=zeros(5*n);M(fids,fids)=MF+MF*B*MS*E*MF;M(fids,cids)=-MF*B*MS;M(cids,fids)=-MS*E*MF;M(cids,cids)=MS;
    rhs=source.rhs;[single,work]=applyRocketSpatialReduction(reduction,coarseFactor,rhs);assert(norm(single(:)-M*rhs(:))/max(1,norm(M*rhs(:)))<1e-9);
    other=reshape(cos((1:5*n)'*.29),n,5);otherValue=applyRocketSpatialReduction(reduction,coarseFactor,other);combined=applyRocketSpatialReduction(reduction,coarseFactor,.7*rhs-1.3*other);
    assert(norm(combined-.7*single+1.3*otherValue,'fro')/max(1,norm(combined,'fro'))<1e-9);
    opt=struct('maximumIterations',5*n,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
    [value,detail]=solveRocketSpatialReductionKrylov(packed,reduction,coarseFactor,rhs,opt);assert(detail.krylov.flag==0&&detail.relativeCResidual<1.1e-10&&norm(value(:)-C\rhs(:))/max(1,norm(C\rhs(:)))<1e-8);
    compact=compactRocketSpatialReductionReport(detail);raw=jsonencode(compact);assert(~contains(raw,'"V":')&&~contains(raw,'"Z":')&&numel(raw)<128*1024);
    [zero,z]=solveRocketSpatialReductionKrylov(packed,reduction,coarseFactor,zeros(n,5),opt);assert(all(zero==0,'all')&&z.CApplications==0&&z.FLocalApplications==0);
    bad=options;bad.maximumReductionBytes=1;rejected=false;try,buildRocketSpatialReduction(packed,source.partition,[],bad);catch ex,rejected=strcmp(ex.identifier,'SpatialReduction:Budget');end;assert(rejected);
    badRepresentatives=split.coarseOriginalCVs;badRepresentatives(1)=badRepresentatives(2);rejected=false;try,buildRocketSpatialSplit(packed,source.partition,badRepresentatives);catch ex,rejected=strcmp(ex.identifier,'SpatialSplit:Representatives');end;assert(rejected);
    cases{j}=struct('packed',packed,'reduction',reduction,'coarseFactor',coarseFactor,'C',C,'A',A,'B',B,'E',E,'D',D,'MF',MF,'S',S,'MS',MS,'M',M,'rhs',rhs,'single',single,'work',work,'value',value,'detail',detail);
end
% 模块三：奇异局部F组只在预条件中退至原CV对角，近似S使用实际替代映射。
toy=struct('kind','complete-five-conservation-block-CSR-unshifted-steady-graph','cellCount',3,'rowPointer',uint32([1;4;7;10]), ...
    'column',uint32(repmat((1:3)',3,1)),'diagonalPosition',uint32([1;5;9]),'order',(1:3)','inverseOrder',(1:3)','values',repmat(eye(5),1,1,9));
part=struct('cellCount',3,'patchCount',1,'members',uint32((1:3)'),'groupID',uint32(ones(3,1)),'sizes',uint8(3),'groupPointer',uint32([1;4]));
fallback=buildRocketSpatialReduction(toy,part,3,options);assert(fallback.Ffactor.factorizationReport.fallbackFBlockCount==1&&norm(densePacked(fallback.schur)+eye(5),'fro')==0);
report=struct('status','passed-two-orderings-spatial-five-field-split-approximate-Schur-back-substitution-and-full-C-reference', ...
    'singularLocalFFallbackKeptExplicit',true,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'spatial_reduction_unit_evidence.mat'),'report','cases','parentFile','toy','fallback','-v7');fid=fopen(fullfile(folder,'spatial_reduction_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function C=densePacked(packed)
n=packed.cellCount;C=zeros(5*n);
for row=1:n,rr=packed.order(row)+(0:4)*n;for edge=double(packed.rowPointer(row)):double(packed.rowPointer(row+1))-1,cc=packed.order(packed.column(edge))+(0:4)*n;C(rr,cc)=packed.values(:,:,edge);end,end
end
function C=denseRect(coupling)
C=zeros(5*coupling.rowCount,5*coupling.columnCount);
for edge=1:numel(coupling.column),rr=double(coupling.owner(edge))+(0:4)*coupling.rowCount;cc=double(coupling.column(edge))+(0:4)*coupling.columnCount;C(rr,cc)=coupling.values(:,:,edge);end
end
function MF=denseFInverse(factor,A)
nf=factor.cellCount;MF=zeros(5*nf);
for k=1:numel(factor.buckets)
    bucket=factor.buckets{k};for page=1:numel(bucket.groupIDs)
        member=double(bucket.members(:,page));ids=reshape(member(:)'+(0:4)'*nf,[],1);block=A(ids,ids);
        if bucket.fallback(page),diagonal=zeros(5*k);for r=1:k,ii=5*(r-1)+(1:5);diagonal(ii,ii)=block(ii,ii);end;block=diagonal;end
        inverse=block\eye(5*k);assert(norm(inverse-bucket.inverse(:,:,page),'fro')/max(1,norm(inverse,'fro'))<1e-8);MF(ids,ids)=inverse;
    end
end
end
function MS=denseCoarseInverse(factor,S)
nc=factor.cellCount;MS=zeros(5*nc);
for k=1:numel(factor.buckets)
    bucket=factor.buckets{k};for page=1:numel(bucket.patchIDs)
        member=double(bucket.members(:,page));ids=reshape(member(:)'+(0:4)'*nc,[],1);block=S(ids,ids);
        if bucket.fallback(page),diagonal=zeros(5*k);for r=1:k,ii=5*(r-1)+(1:5);diagonal(ii,ii)=block(ii,ii);end;block=diagonal;end
        MS(ids,ids)=block\eye(5*k);
    end
end
end
