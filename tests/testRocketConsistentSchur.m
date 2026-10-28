function report=testRocketConsistentSchur(folder)
%TESTROCKETCONSISTENTSCHUR 对两种排序核对完整固定动量逆Schur与稠密参考。
% 模块一：从既有完整小域C及动量因子建立作用，独立核对D-E*(LU)^{-1}*B。
% 模块二：保存各中间量以便Python复算，检查线性及非法输入，不增加物理解面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'ConsistentSchurTest:Overwrite','不得覆盖一致Schur测试。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','fieldschurunit01','field_schur_unit_evidence.mat');p=load(parentFile,'cases');cases=cell(2,1);
for j=1:2
    c=p.cases{j};a=load(fullfile(root,'outputs','v83_validation','schurcheckpointunit01',sprintf('assembly%d',j),'schur_assembly_evidence.mat'));
    m=load(fullfile(root,'outputs','v83_validation','schurcheckpointunit01',sprintf('assembly%d',j),'momentum_factor_evidence.mat'),'M');n=m.M.cellCount;
    packed=c.packed;operator=struct('kind','fixed-momentum-inverse-consistent-Schur-action','cellCount',n, ...
        'rowPointer',packed.rowPointer,'column',packed.column,'owners',uint32(repelem((1:n)',diff(double(packed.rowPointer)))), ...
        'D',packed.values([1 5],[1 5],:),'B',a.couplingB,'E',a.couplingE,'momentumFactor',m.M);
    value=reshape(sin((1:2*n)'*.37),n,2);[image,detail]=applyRocketConsistentSchur(operator,value);
    [observed,observation]=observeRocketFieldILU0(m.M,detail.Bvalue);assert(observation.valueIsUsable&&isequaln(observed,detail.momentumValue));
    exactApproximateMomentumInverse=c.MU\(c.ML\eye(3*n));reference=c.D-c.E*exactApproximateMomentumInverse*c.B;
    relative=norm(image(:)-reference*value(:))/max(1,norm(reference*value(:)));assert(relative<1e-10);
    assert(norm(detail.momentumValue(:)-exactApproximateMomentumInverse*c.B*value(:))/max(1,norm(detail.momentumValue,'fro'))<1e-10);
    other=reshape(cos((1:2*n)'*.23),n,2);otherImage=applyRocketConsistentSchur(operator,other);combined=applyRocketConsistentSchur(operator,.7*value-1.3*other);
    linearity=norm(combined-.7*image+1.3*otherImage,'fro')/max(1,norm(combined,'fro'));assert(linearity<1e-10);
    rejected=false;try,applyRocketConsistentSchur(operator,zeros(n,3));catch ex,rejected=strcmp(ex.identifier,'ConsistentSchur:Input');end;assert(rejected);
    cases{j}=struct('operator',operator,'value',value,'image',image,'detail',detail,'reference',reference,'relativeDifference',relative,'linearity',linearity,'observation',observation);
end
% 模块三：合成上下三角溢出及对角乘溢出，旧入口仍拒绝，观察器输出不得伪装成零。
toy=struct('kind','factored-conservative-field-block-ILU0','dimension',2,'cellCount',3,'rowPointer',uint32([1;2;4;6]), ...
    'column',uint32([1;1;2;2;3]),'diagonalPosition',uint32([1;3;5]),'values',repmat(eye(2),1,1,5),'inverseDiagonal',repmat(eye(2),1,1,3));
toy.values(:,:,[2 4])=repmat(1e200*eye(2),1,1,2);[bad,lower]=observeRocketFieldILU0(toy,ones(3,2));assert(~lower.valueIsUsable&&all(isnan(bad),'all')&&lower.firstFailure.algebraicRow==3&&strcmp(lower.firstFailure.stage,'forward-lower-sum'));
rejected=false;try,solveRocketFieldILU0(toy,ones(3,2));catch ex,rejected=strcmp(ex.identifier,'FieldILU:Nonfinite');end;assert(rejected);
upperToy=toy;upperToy.rowPointer=uint32([1;3;5;6]);upperToy.column=uint32([1;2;2;3;3]);upperToy.diagonalPosition=uint32([1;3;5]);upperToy.values=repmat(eye(2),1,1,5);upperToy.values(:,:,[2 4])=repmat(1e200*eye(2),1,1,2);
[bad,upper]=observeRocketFieldILU0(upperToy,ones(3,2));assert(~upper.valueIsUsable&&all(isnan(bad),'all')&&upper.firstFailure.algebraicRow==1&&strcmp(upper.firstFailure.stage,'backward-upper-sum'));
diagonalToy=toy;diagonalToy.cellCount=1;diagonalToy.rowPointer=uint32([1;2]);diagonalToy.column=uint32(1);diagonalToy.diagonalPosition=uint32(1);diagonalToy.values=1e-200*eye(2);diagonalToy.inverseDiagonal=1e200*eye(2);
[bad,diagonalFailure]=observeRocketFieldILU0(diagonalToy,1e200*ones(1,2));assert(~diagonalFailure.valueIsUsable&&all(isnan(bad),'all')&&strcmp(diagonalFailure.firstFailure.stage,'backward-diagonal-multiply'));
% 模块四：独立SVD核对近等奇异值、近秩亏、零块和极宽尺度块的诊断。
diagonalTests=cat(3,eye(2),[1 1e-8;-1e-8 1+1e-10],[1 1;1 1+1e-10],zeros(2),1e250*[1 .2;-.3 .9],1e-250*[1 .2;-.3 .9]);
ratios=rocketTwoFieldSingularRatio(diagonalTests);referenceRatios=zeros(size(ratios));for k=1:numel(ratios),sv=svd(diagonalTests(:,:,k));if sv(1)>0,referenceRatios(k)=sv(2)/sv(1);end,end;assert(max(abs(ratios-referenceRatios))<1e-14);
report=struct('status','passed-small-consistent-fixed-momentum-inverse-Schur-linearity-and-first-overflow-observation','cases',2, ...
    'finiteObserverExactlyMatchesDefault',true,'threeFirstOverflowStagesVerified',true,'oldNonfiniteRejectionPreserved',true, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'consistent_schur_unit_evidence.mat'),'report','cases','parentFile','toy','lower','upperToy','upper','diagonalToy','diagonalFailure','diagonalTests','ratios','referenceRatios','-v7');
fid=fopen(fullfile(folder,'consistent_schur_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
