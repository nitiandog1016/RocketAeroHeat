function report=testRocketPeakWeightedDirection(folder)
%TESTROCKETPEAKWEIGHTEDDIRECTION 独立稠密检验冻结峰值权重、完整转置和增广总方向。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);s=load(fullfile(root,'outputs','v83_validation','updateradiusunit01','regularized_parent_unit_evidence.mat'),'cases');
p=load(fullfile(root,'outputs','v83_validation','parentlsqrunit02','parent_least_squares_unit_evidence.mat'),'cases');cases=cell(2,1);
for j=1:2
    a=s.cases{j};b=p.cases{j};n=size(a.Q,1);N=5*n;
    reference=reshape(20*sin((1:N)'),n,5);[metric,weighting]=reweightRocketPrimitiveUpdateMetric(a.metric,reference,struct('transitionFraction',.2,'maximumRowWeight',64));
    expectedWeights=min(64,sqrt(1+(reshape(a.W*reference(:),n,5)/.2).^2));
    assert(max(abs(metric.peakRowWeights-expectedWeights),[],'all')<1e-12);
    W=expectedWeights(:).*a.W;Wactual=zeros(N);WTactual=Wactual;
    for k=1:N
        v=zeros(n,5);v(k)=1;z=applyRocketPrimitiveUpdateMetric(metric,v);Wactual(:,k)=z(:);
        z=applyRocketPrimitiveUpdateMetric(metric,v,'transp');WTactual(:,k)=z(:);
    end
    assert(norm(W-Wactual,'fro')<1e-12*norm(W,'fro')&&norm(W'-WTactual,'fro')<1e-12*norm(W,'fro'));
    zero=reweightRocketPrimitiveUpdateMetric(a.metric,zeros(n,5));v=reshape(cos((1:N)'),n,5);
    assert(isequal(applyRocketPrimitiveUpdateMetric(zero,v),applyRocketPrimitiveUpdateMetric(a.metric,v)));
    [value,detail,coordinates]=solveRocketRegularizedParentDirection(b.packed,b.factor,a.rhs,a.initial,a.columnScale,metric, ...
        struct('maximumIterations',900,'regularizationStrength',.1));
    A=[a.C;detail.penaltyInOriginalResidualUnits*W];dense=A\[a.rhs(:);zeros(N,1)];
    difference=norm(value(:)-dense)/max(1,norm(dense));assert(difference<1e-8&&detail.adjointIdentityRelativeDifference<1e-12);
    cases{j}=struct('base',a,'reference',reference,'metric',metric,'weighting',weighting,'W',W,'value',value, ...
        'detail',detail,'coordinates',coordinates,'denseDifference',difference);
end
report=struct('status','passed-two-orderings-frozen-peak-rows-full-adjoint-and-dense-total-direction', ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'peak_weighted_unit_evidence.mat'),'report','cases','-v7');fprintf('%s\n',jsonencode(report));
end
