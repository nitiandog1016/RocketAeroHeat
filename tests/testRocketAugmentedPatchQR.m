function report=testRocketAugmentedPatchQR(folder)
%TESTROCKETAUGMENTEDPATCHQR 用两种排列、完整稠密增广解和全部局部Gram校核新右映射。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
p=load(fullfile(root,'outputs/v83_validation/parentlsqrunit02/parent_least_squares_unit_evidence.mat'),'cases');
m=load(fullfile(root,'outputs/v83_validation/peakweightunit01/peak_weighted_unit_evidence.mat'),'cases');cases=cell(2,1);
for j=1:2
    a=p.cases{j};b=m.cases{j};n=a.packed.cellCount;N=5*n;penalty=b.detail.penaltyInOriginalResidualUnits;
    factor=buildRocketAugmentedPatchQR(a.packed,a.factor,b.metric,penalty);P=zeros(N);PT=P;
    for column=1:N
        e=zeros(n,5);e(column)=1;v=applyRocketAugmentedPatchQR(factor,e);P(:,column)=v(:);
        v=applyRocketAugmentedPatchQR(factor,e,'transp');PT(:,column)=v(:);
    end
    assert(norm(PT-P','fro')<1e-12*norm(P,'fro'));
    rhs=b.base.rhs;dense=[b.base.C;penalty*b.W]\[rhs(:);zeros(N,1)];
    [value,detail,coordinates]=solveRocketAugmentedQRDirection(a.packed,factor,rhs,b.metric,struct('maximumIterations',900));
    error=norm(value(:)-dense)/max(1,norm(dense));assert(error<1e-8&&detail.adjointIdentityRelativeDifference<1e-12);
    assert(norm(P*coordinates(:)-value(:))<1e-12*max(1,norm(value(:))));
    [warm,warmReport]=solveRocketAugmentedQRDirection(a.packed,factor,rhs,b.metric,struct('maximumIterations',900),b.base.initial);
    warmError=norm(warm(:)-dense)/max(1,norm(dense));assert(warmError<1e-8);
    cases{j}=struct('source',b,'packed',a.packed,'factor',factor,'P',P,'rhs',rhs,'dense',dense, ...
        'value',value,'detail',detail,'coordinates',coordinates,'warm',warm,'warmReport',warmReport, ...
        'denseDifference',error,'warmDenseDifference',warmError);
end
report=struct('status','passed-two-orderings-all-map-columns-transpose-zero-and-warm-dense-total-objective', ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'physicalSolution',false);
save(fullfile(folder,'augmented_qr_unit_evidence.mat'),'cases','report','-v7');disp(report);
end
