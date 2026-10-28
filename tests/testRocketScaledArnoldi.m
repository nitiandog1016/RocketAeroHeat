function report=testRocketScaledArnoldi(folder)
%TESTROCKETSCALEDARNOLDI 验证显式整方向缩放、物理像、旧默认与原拒绝门槛。
% 模块一：旧完整增强内核逐值回放；含失败观察器的非线性拒绝仍生效。
% 模块二：跨32数量级的可变标量右预条件，独立稠密作用和最小剩余验证。
% 模块三：小域完整五量C使用已有因子，禁止新增当前J/NS或流场接受声明。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'ScaledArnoldiTest:Overwrite','不得覆盖测试证据。');mkdir(folder);
testRocketKrylovFailureEvidence(fullfile(folder,'default_replay'));
n=48;A=diag(linspace(1,8,n))+diag(.4*ones(n-1,1),1)-diag(.2*ones(n-1,1),-1);
b=sin((1:n)'*.27);gains=10.^linspace(-16,16,n)';calls=0;
options=struct('restartLength',48,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',48, ...
    'absoluteTolerance',1e-10,'maximumBasisBytes',32*1024^2,'saveBasisEvidence',true,'robustArnoldi',true);
[value,detail]=solveRocketRecycledGMRES(@(v)A*v,@gain,b,options);direct=A\b;
assert(detail.flag==0&&norm(A*value-b)/norm(b)<1e-9&&norm(value-direct)/norm(direct)<1e-9);
assert(detail.maximumBasisOrthogonalityError<1e-8&&detail.originalModelRelationThreshold==1e-8);
entry=detail.cycles{1};used=entry.arnoldiSteps;expected=zeros(n,used);
for j=1:used,expected(:,j)=gains(j)*entry.V(:,j)/entry.imageNormalizationScales(j);end
assert(norm(entry.Z-expected,'fro')<1e-12&&norm(A*entry.Z-entry.V*entry.smallH,'fro')<1e-11);
assert(max(entry.imageNormalizationScales)/min(entry.imageNormalizationScales)>1e12);
zeroOptions=options;zeroOptions.restartLength=4;zeroOptions.maximumArnoldiSteps=4;
[nullValue,nullReport]=solveRocketRecycledGMRES(@(v)zeros(size(v)),@(v)v,ones(8,1),zeroOptions);
assert(all(nullValue==0)&&nullReport.flag~=0&&nullReport.cycles{1}.imageNormalizationScales==1&& ...
    strcmp(nullReport.stopReason,'Arnoldi-breakdown-with-unresolved-current-model-residual'));
s=load(fullfile(root,'outputs','v83_validation','blockiluunit01','block_ilu_unit_evidence.mat'));
operator=s.operator;factor=s.cases{2}.factor;C=s.C;rhs=s.cases{2}.rhs;
innerOptions=struct('maximumIterations',64,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true,'robustArnoldi',true);
[smallValue,smallReport]=solveRocketBlockGraphKrylov(operator,factor,rhs,innerOptions);
assert(smallReport.krylov.flag==0&&norm(C*smallValue(:)-rhs(:))/norm(rhs(:))<1.1e-10);
report=struct('status','passed-explicit-scaled-adaptive-Arnoldi-default-replay-variable-gain-and-small-full-C', ...
    'variableGainArnoldiSteps',detail.arnoldiSteps,'relativeResidual',norm(A*value-b)/norm(b), ...
    'basisOrthogonalityError',detail.maximumBasisOrthogonalityError,'smallRelativeCResidual',smallReport.relativeCResidual,'zeroImageRejectedAsUnresolved',true, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'scaled_arnoldi_unit_evidence.mat'),'report','A','b','gains','value','detail','direct','C','rhs','smallValue','smallReport','factor','nullValue','nullReport','-v7');
fid=fopen(fullfile(folder,'scaled_arnoldi_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function z=gain(v)
        calls=calls+1;z=gains(calls)*v;
    end
end
