function report=testRocketBlockFillLevel1(folder)
%TESTROCKETBLOCKFILLLEVEL1 小域两排序的一级填充、原矩阵与完整三角解验证。
% 模块一：独立稠密图乘积核对全部符号位置，原块逐值及新增零块逐值检查。
% 模块二：扩充模式上的零填充分解实现一级块ILU，完整L/U可独立复算。
% 模块三：预算与二次填充拒绝；不调用新当前J/原NS。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'BlockFillTest:Overwrite','不得覆盖一级填充测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','blockiluunit01','block_ilu_unit_evidence.mat');s=load(sourceFile);C=s.C;cases=cell(2,1);errors=zeros(2,2);
for j=1:2
    original=s.cases{j}.packed;before=original;enriched=enrichRocketBlockGraphLevel1(original);e=enriched.enrichmentReport;n=original.cellCount;
    P=zeros(n);for row=1:n,ids=double(original.rowPointer(row)):double(original.rowPointer(row+1))-1;P(row,original.column(ids))=1;end
    expected=(P+tril(P,-1)*triu(P,1))>0;actual=false(n);
    for row=1:n,ids=double(enriched.rowPointer(row)):double(enriched.rowPointer(row+1))-1;actual(row,enriched.column(ids))=true;end
    assert(isequal(expected,actual)&&isequaln(original,before));
    assert(isequaln(original.values,enriched.values(:,:,e.originalToEnrichedPosition)));
    added=true(1,numel(enriched.column));added(e.originalToEnrichedPosition)=false;assert(all(enriched.values(:,:,added)==0,'all'));
    probe=s.cases{j}.probe;value=applyRocketPackedGraph(enriched,probe);errors(j,1)=norm(value(:)-C*probe(:))/max(1,norm(C*probe(:)));
    twice=false;try,enrichRocketBlockGraphLevel1(enriched);catch ex,twice=strcmp(ex.identifier,'BlockFill:Kind');end;assert(twice);
    factor=factorRocketBlockILU0(enriched);[L,U]=matrices(factor);rhs=s.cases{j}.rhs;solution=solveRocketBlockILU0(factor,rhs);
    ordered=rhs(factor.order,:)';direct=U\(L\ordered(:));observed=solution(factor.order,:)';errors(j,2)=norm(observed(:)-direct)/max(1,norm(direct));
    assert(max(errors(j,:))<1e-9&&e.persistentFactorStorageUpperBoundBytes<=original.options.maximumPersistentBytes);
    cases{j}=struct('original',original,'enriched',enriched,'factor',factor,'L',L,'U',U,'rhs',rhs,'solution',solution);
end
bad=original;bad.options.maximumPersistentBytes=1;rejected=false;try,enrichRocketBlockGraphLevel1(bad);catch ex,rejected=strcmp(ex.identifier,'BlockFill:Budget');end;assert(rejected);
report=struct('status','passed-small-two-orderings-symbolic-level1-original-values-zero-fill-and-dense-triangular-reference', ...
    'relativeErrors',errors,'budgetRejected',rejected,'doubleEnrichmentRejected',twice,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'block_fill1_unit_evidence.mat'),'report','sourceFile','C','cases','-v7');
fid=fopen(fullfile(folder,'block_fill1_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function [L,U]=matrices(factor)
n=factor.cellCount;L=eye(5*n);U=zeros(5*n);
for row=1:n
    ids=double(factor.rowPointer(row)):double(factor.rowPointer(row+1))-1;rr=(row-1)*5+(1:5);
    for entry=ids
        col=double(factor.column(entry));cc=(col-1)*5+(1:5);
        if col<row,L(rr,cc)=factor.values(:,:,entry);else,U(rr,cc)=factor.values(:,:,entry);end
    end
end
end
