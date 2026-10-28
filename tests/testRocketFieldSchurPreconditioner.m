function report=testRocketFieldSchurPreconditioner(folder)
%TESTROCKETFIELDSCHURPRECONDITIONER 两排序完整小域Schur装配/块三角解/固定线性核对。
% 模块一：冻结完整C独立按字段构造S=D-E*diag(A)^{-1}*B，全部原CV均参与。
% 模块二：上三角的完整2/3量对角块均由稠密解核对，不作标量三角假设。
% 模块三：原C保持，预算/非法右端拒绝，无新J/NS或覆盖。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FieldSchurTest:Overwrite','不得覆盖Schur单位证据。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','blockiluunit01','block_ilu_unit_evidence.mat');s=load(sourceFile);C=s.C;n=size(s.Q,1);cases=cell(2,1);errors=zeros(2,3);
for j=1:2
    packed=s.cases{j}.packed;before=packed;pre=buildRocketFieldSchurPreconditioner(packed);assert(isequaln(packed,before));order=packed.order;
    momentum=reshape(order+n*(1:3),[],1);scalar=reshape(order+n*[0 4],[],1);A=C(momentum,momentum);B=C(momentum,scalar);E=C(scalar,momentum);D=C(scalar,scalar);local=zeros(3*n);
    for row=1:n,ids=row+n*(0:2);local(ids,ids)=A(ids,ids)\eye(3);end
    referenceS=D-E*local*B;[ML,MU]=matrices(pre.momentumFactor);[SL,SU]=matrices(pre.schurFactor);
    [~,mmask]=matrixShape(pre.momentumFactor);[~,smask]=matrixShape(pre.schurFactor);
    % matrices输出按字段分组，以便直接与C的上述子块对应。
    errors(j,1)=norm((ML*MU-A).*mmask,'fro')/max(1,norm(A,'fro'));
    errors(j,2)=norm((SL*SU-referenceS).*smask,'fro')/max(1,norm(referenceS,'fro'));
    rhs=s.cases{j}.rhs;[value,detail]=solveRocketFieldSchurPreconditioner(pre,rhs);r=rhs(order,:);rm=reshape(r(:,2:4),[],1);rs=reshape(r(:,[1 5]),[],1);
    ym=MU\(ML\rm);ys=SU\(SL\(rs-E*ym));ym=ym-MU\(ML\(B*ys));expected=zeros(n,5);expected(:,2:4)=reshape(ym,n,3);expected(:,[1 5])=reshape(ys,n,2);physical=zeros(n,5);physical(order,:)=expected;
    errors(j,3)=norm(value-physical,'fro')/max(1,norm(physical,'fro'));assert(max(errors(j,:))<1e-8,'FieldSchurTest:Dense','完整Schur装配/前后代未通过稠密参考。');
    other=reshape(cos((1:numel(rhs))'*.29),size(rhs));a=solveRocketFieldSchurPreconditioner(pre,other);combined=solveRocketFieldSchurPreconditioner(pre,.7*rhs-1.3*other);
    linearity=norm(combined-.7*value+1.3*a,'fro')/max(1,norm(combined,'fro'));assert(linearity<1e-9);
    cases{j}=struct('pre',pre,'packed',packed,'A',A,'B',B,'E',E,'D',D,'referenceS',referenceS,'ML',ML,'MU',MU,'SL',SL,'SU',SU,'rhs',rhs,'value',value,'linearity',linearity,'detail',detail);
end
bad=packed;bad.options.maximumPersistentBytes=1;rejected=false;try,buildRocketFieldSchurPreconditioner(bad);catch ex,rejected=strcmp(ex.identifier,'FieldSchur:Budget');end;assert(rejected);
badRhs=false;try,solveRocketFieldSchurPreconditioner(pre,zeros(n,4));catch ex,badRhs=strcmp(ex.identifier,'FieldSchur:RHS');end;assert(badRhs);
report=struct('status','passed-small-two-orderings-complete-Schur-field-assembly-triangular-linearity-and-budget-reference', ...
    'relativeErrors',errors,'budgetRejected',rejected,'badRhsRejected',badRhs,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'field_schur_unit_evidence.mat'),'report','sourceFile','C','cases','-v7');
fid=fopen(fullfile(folder,'field_schur_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function [L,U]=matrices(f)
n=f.cellCount;d=f.dimension;L=eye(n*d);U=zeros(n*d);
for row=1:n
    for entry=double(f.rowPointer(row)):double(f.rowPointer(row+1))-1
        col=double(f.column(entry));rr=row+n*(0:d-1);cc=col+n*(0:d-1);
        if col<row,L(rr,cc)=f.values(:,:,entry);else,U(rr,cc)=f.values(:,:,entry);end
    end
end
end
function [dummy,mask]=matrixShape(f)
n=f.cellCount;d=f.dimension;mask=zeros(n*d);dummy=[];
for row=1:n
    ids=double(f.rowPointer(row)):double(f.rowPointer(row+1))-1;
    for col=double(f.column(ids))',mask(row+n*(0:d-1),col+n*(0:d-1))=1;end
end
end
