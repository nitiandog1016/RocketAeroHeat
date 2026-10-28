function report=testRocketFullFImageScale(folder)
%TESTROCKETFULLFIMAGESCALE 核对全部完整像列、排列、缩放伴随、LSQR及旧缺省数值路径。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);sourceFile=fullfile(root,'outputs','v83_validation','fullFlsqrunit02','full_F_least_squares_unit_evidence.mat');s=load(sourceFile,'cases');cases=cell(2,1);
for j=1:2
    a=s.cases{j};graph=a.graph;factor=a.factor;
    [scale,norms,construction]=buildRocketFullFImageScale(graph,factor);expected=reshape(vecnorm(a.A*a.M,2,1),size(norms));
    assert(max(abs(norms-expected)./expected,[],'all')<1e-12);
    [permutedScale,permutedNorms]=buildRocketFullFImageScale(a.permutedGraph,factor);permutedExpected=reshape(vecnorm(a.permutedA*a.M,2,1),size(norms));
    assert(max(abs(permutedNorms-permutedExpected)./permutedExpected,[],'all')<1e-12);
    rhs=a.rhs;initial=a.initial;options=struct('maximumIterations',300,'relativeTolerance',1e-10);
    [value,detail]=solveRocketFullFLeastSquares(graph,factor,rhs,initial,options,scale);direct=a.A\rhs(:);
    assert(detail.fullFInternalTargetMet&&detail.columnScalingUsed&&detail.adjointIdentityRelativeDifference<1e-12&&norm(value(:)-direct)/max(1,norm(direct))<1e-8);
    [plain,plainDetail]=solveRocketFullFLeastSquares(graph,factor,rhs,initial,options);
    [identity,identityDetail]=solveRocketFullFLeastSquares(graph,factor,rhs,initial,options,ones(size(rhs)));
    assert(isequal(plain,identity)&&isequal(plain,a.value));
    assert(isequal(plainDetail.nativeResidualHistory,a.detail.nativeResidualHistory)&&isequal(plainDetail.operatorApplications,a.detail.operatorApplications));
    assert(isequal(plainDetail.nativeResidualHistory,identityDetail.nativeResidualHistory));
    rejected=false;try,buildRocketFullFImageScale(graph,factor,struct('maximumWorkspaceBytes',1));catch ex,rejected=strcmp(ex.identifier,'FImageScale:Budget');end;assert(rejected);
    bad=scale;bad(1)=0;rejected=false;try,solveRocketFullFLeastSquares(graph,factor,rhs,initial,options,bad);catch ex,rejected=strcmp(ex.identifier,'FullFLSQR:Scale');end;assert(rejected);
    cases{j}=struct('graph',graph,'factor',factor,'A',a.A,'M',a.M,'columnScale',scale,'columnNorm',norms,'construction',construction,'rhs',rhs,'initial',initial,'value',value,'detail',detail, ...
        'permutedGraph',a.permutedGraph,'permutedA',a.permutedA,'permutedScale',permutedScale,'permutedNorms',permutedNorms,'plainValue',plain,'plainDetail',plainDetail);
end
report=struct('status','passed-all-small-image-columns-permutations-scaled-adjoint-LSQR-and-old-default-replay', ...
    'cases',2,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'full_F_image_scale_unit_evidence.mat'),'report','cases','sourceFile','-v7');fprintf('%s\n',jsonencode(report));
end
