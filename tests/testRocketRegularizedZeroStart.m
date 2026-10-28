function report=testRocketRegularizedZeroStart(folder)
%TESTROCKETREGULARIZEDZEROSTART 固定同一增广目标，比较非零与零初猜的独立稠密解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);s=load(fullfile(root,'outputs','v83_validation','peakweightunit01','peak_weighted_unit_evidence.mat'),'cases');
p=load(fullfile(root,'outputs','v83_validation','parentlsqrunit02','parent_least_squares_unit_evidence.mat'),'cases');cases=cell(2,1);
for j=1:2
    a=s.cases{j};b=p.cases{j};n=size(a.base.Q,1);options=struct('maximumIterations',900,'regularizationStrength',.1,'metricNormalization',a.detail.metricNormalization);
    [zero,zeroReport,zeroCoordinates]=solveRocketRegularizedParentDirection(b.packed,b.factor,a.base.rhs,zeros(n,5),a.base.columnScale,a.metric,options);
    [warm,warmReport]=solveRocketRegularizedParentDirection(b.packed,b.factor,a.base.rhs,a.base.initial,a.base.columnScale,a.metric,options);
    assert(zeroReport.penaltyInOriginalResidualUnits==warmReport.penaltyInOriginalResidualUnits&& ...
        zeroReport.metricNormalization==warmReport.metricNormalization&&zeroReport.zeroDirectionNormalizedObjective==1);
    A=[a.base.C;zeroReport.penaltyInOriginalResidualUnits*a.W];dense=A\[a.base.rhs(:);zeros(5*n,1)];
    assert(norm(zero(:)-dense)/max(1,norm(dense))<1e-8&&norm(warm(:)-dense)/max(1,norm(dense))<1e-8);
    raw=applyRocketPrimitiveUpdateMetric(a.base.metric,zero);
    assert(abs(zeroReport.maximumDensityFraction-max(abs(raw(:,1))))<1e-12*max(1,zeroReport.maximumDensityFraction));
    assert(zeroReport.primitiveUpdateNormIncludesRowWeights&&zeroReport.beatsZeroDirectionObjective);
    cases{j}=struct('source',a,'zeroValue',zero,'warmValue',warm,'zeroReport',zeroReport,'warmReport',warmReport, ...
        'zeroCoordinates',zeroCoordinates,'dense',dense,'unweightedFinal',raw);
end
report=struct('status','passed-fixed-augmented-objective-zero-and-warm-dense-solutions-and-raw-state-labels', ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'peak_zero_unit_evidence.mat'),'cases','report','-v7');fprintf('%s\n',jsonencode(report));
end
