function report=testRocketPatchHierarchy(folder)
%TESTROCKETPATCHHIERARCHY 两排序核对全部系数RCP、R/P和固定V循环稠密参考。
% 模块一：明确代数平均限制不等于体积守恒模型，原C与物理CV不变。
% 模块二：完整V循环、线性、原C求解及全部物理Arnoldi方向共同验收。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'PatchHierarchyTest:Overwrite','不得覆盖完整层次测试。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','coupledpatchunit03','coupled_patches_unit_evidence.mat');s=load(parentFile,'cases','toy','toyFactor');cases=cell(2,1);
for j=1:2
    source=s.cases{j};packed=source.packed;factor=source.factor;n=packed.cellCount;C=denseGraph(packed);
    options=struct('maximumHierarchyBytes',32*1024^2,'maximumConstructionBytes',64*1024^2);
    hierarchy=buildRocketPatchHierarchy(packed,factor,options);assert(isequaln(packed,source.packed)&&isequaln(factor,source.factor));
    assert(numel(hierarchy.levels)==2&&hierarchy.levels{2}.terminal);
    partition=factor.partition;ng=partition.patchCount;P=zeros(n,ng);P(sub2ind(size(P),(1:n)',double(partition.groupID)))=1;
    R=P'./double(partition.sizes(:));P5=kron(eye(5),P);R5=kron(eye(5),R);coarse=denseGraph(hierarchy.levels{2}.packed);
    assert(norm(R5*P5-eye(5*ng),'fro')<1e-13&&norm(coarse-R5*C*P5,'fro')/max(1,norm(coarse,'fro'))<1e-12);
    M=source.P+P5*(coarse\(R5*(eye(5*n)-C*source.P)));
    rhs=source.rhs;[single,application]=applyRocketPatchHierarchy(hierarchy,rhs);assert(norm(single(:)-M*rhs(:))/max(1,norm(M*rhs(:)))<1e-9);
    other=reshape(cos((1:5*n)'*.29),n,5);otherValue=applyRocketPatchHierarchy(hierarchy,other);combined=applyRocketPatchHierarchy(hierarchy,.7*rhs-1.3*other);
    assert(norm(combined-.7*single+1.3*otherValue,'fro')/max(1,norm(combined,'fro'))<1e-9);
    opt=struct('maximumIterations',5*n,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
    [value,detail]=solveRocketPatchHierarchyKrylov(hierarchy,rhs,opt);assert(detail.krylov.flag==0&&detail.relativeCResidual<1.1e-10);
    direct=C\rhs(:);assert(norm(value(:)-direct)/max(1,norm(direct))<1e-8);
    compact=compactRocketPatchHierarchyReport(detail);raw=jsonencode(compact);assert(~contains(raw,'"V":')&&~contains(raw,'"Z":')&&numel(raw)<128*1024);
    assert(detail.preconditionerCApplicationsByLevel(1)==detail.VcycleApplications&&detail.preconditionerTerminalApplicationsByLevel(2)==detail.VcycleApplications);
    assert(detail.totalRootCApplications==detail.solverCApplications+detail.VcycleApplications);
    [zero,z]=solveRocketPatchHierarchyKrylov(hierarchy,zeros(n,5),opt);assert(all(zero==0,'all')&&z.totalRootCApplications==0&&z.VcycleApplications==0);
    bad=options;bad.maximumHierarchyBytes=1;rejected=false;try,buildRocketPatchHierarchy(packed,factor,bad);catch ex,rejected=strcmp(ex.identifier,'PatchHierarchy:Budget');end;assert(rejected);
    badK=opt;badK.maximumBasisBytes=1;rejected=false;try,solveRocketPatchHierarchyKrylov(hierarchy,rhs,badK);catch ex,rejected=strcmp(ex.identifier,'RecycledGMRES:Budget');end;assert(rejected);
    deepOptions=options;deepOptions.terminalMaximumCells=1;deepHierarchy=buildRocketPatchHierarchy(packed,factor,deepOptions);
    assert(numel(deepHierarchy.levels)>=3);deepMatrix=denseHierarchyReference(deepHierarchy);[deepSingle,deepWork]=applyRocketPatchHierarchy(deepHierarchy,rhs);
    assert(norm(deepSingle(:)-deepMatrix*rhs(:))/max(1,norm(deepMatrix*rhs(:)))<1e-9);
    assert(all(deepWork.CApplicationsByLevel(1:end-1)==1)&&deepWork.CApplicationsByLevel(end)==0&&deepWork.terminalApplicationsByLevel(end)==1);
    cases{j}=struct('packed',packed,'hierarchy',hierarchy,'C',C,'P5',P5,'R5',R5,'coarse',coarse,'M',M,'rhs',rhs,'single',single,'application',application,'value',value,'detail',detail, ...
        'deepHierarchy',deepHierarchy,'deepMatrix',deepMatrix,'deepSingle',deepSingle,'deepWork',deepWork);
end
% 模块四：最粗完整系统奇异时明确拒绝，不以原CV对角替代终端系统。
singularRejected=false;
try,buildRocketPatchHierarchy(s.toy,s.toyFactor,struct('failureEvidenceFile',fullfile(folder,'singular_terminal_evidence.mat')));
catch ex,singularRejected=strcmp(ex.identifier,'PatchHierarchy:TerminalInverse');end
assert(singularRejected);
report=struct('status','passed-two-orderings-all-coefficient-RCP-fixed-V-cycle-and-complete-C-reference', ...
    'deepRecursiveCycleTested',true,'singularTerminalRejected',singularRejected,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'patch_hierarchy_unit_evidence.mat'),'report','cases','parentFile','-v7');fid=fopen(fullfile(folder,'patch_hierarchy_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function M=denseHierarchyReference(hierarchy)
% 模块五：逐层稠密乘法独立组成完整固定V循环，不调用递归作用入口。
for level=numel(hierarchy.levels):-1:1
    entry=hierarchy.levels{level};C=denseGraph(entry.packed);n=entry.packed.cellCount;
    if entry.terminal,M=C\eye(5*n);continue;end
    part=entry.factor.partition;ng=part.patchCount;P=zeros(n,ng);P(sub2ind(size(P),(1:n)',double(part.groupID)))=1;
    R=P'./double(part.sizes(:));P5=kron(eye(5),P);R5=kron(eye(5),R);local=zeros(5*n);
    for k=1:numel(entry.factor.buckets)
        bucket=entry.factor.buckets{k};
        for page=1:numel(bucket.patchIDs)
            member=double(bucket.members(:,page));ids=reshape(member(:)'+(0:4)'*n,[],1);block=C(ids,ids);
            if bucket.fallback(page),error('PatchHierarchyTest:UnexpectedFallback','此非奇异小域不应退至CV对角。');end
            local(ids,ids)=block\eye(5*k);
        end
    end
    M=local+P5*M*R5*(eye(5*n)-C*local);
end
end
function C=denseGraph(packed)
% 模块三：独立逐块原物理编号稠密装配，不调用聚合构造器。
n=packed.cellCount;C=zeros(5*n);
for row=1:n
    rr=packed.order(row)+(0:4)*n;
    for edge=double(packed.rowPointer(row)):double(packed.rowPointer(row+1))-1
        cc=packed.order(packed.column(edge))+(0:4)*n;C(rr,cc)=packed.values(:,:,edge);
    end
end
end
