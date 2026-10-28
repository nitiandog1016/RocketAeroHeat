function [value,report,coordinates]=solveRocketAugmentedQRDirection(packed,factor,rhs,metric,options,initial)
%SOLVEROCKETAUGMENTEDQRDIRECTION 在同一完整增广目标上用QR右预处理求方向。
% 模块一：min ||C*v-rhs||²+penalty²||W*v||²；完整五量及原残差范数保持。
% 模块二：LSQR仅改变v=initial+P*y中的P；固定penalty，非零初猜惩罚总方向。
% 模块三：分别报告原C1e-3、优化目标、真实状态幅度；不授予原NS或物理资格。
if nargin<5,options=struct;end
if nargin<6,initial=zeros(size(rhs));end
defaults=struct('maximumIterations',128,'relativeTolerance',1e-3);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'AugmentedQR:Options','未知QR方向求解选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(options.relativeTolerance,{'numeric'},{'scalar','real','finite','positive','<=',1e-3});
shape=[packed.cellCount 5];N=prod(shape);penalty=factor.penaltyInOriginalResidualUnits;
assert(isequal(size(rhs),shape)&&isequal(size(initial),shape)&&metric.cellCount==packed.cellCount&& ...
    all(isfinite(rhs),'all')&&all(isfinite(initial),'all')&&8*N*32<=384*1024^2, ...
    'AugmentedQR:Input','完整方向输入非法或额外向量数组超384MiB声明预算。');
timer=tic;counts=struct('C',0,'CT',0,'P',0,'PT',0,'W',0,'WT',0);
denominator=max(realmin,norm(rhs,'fro'));initialImage=original(initial);initialMetric=weight(initial,'notransp');
augmentedRhs=[reshape(rhs-initialImage,[],1);-penalty*initialMetric(:)];
x=sin((1:N)');z=cos((1:2*N)'/sqrt(2));bx=operate(x,'notransp');btz=operate(z,'transp');
adjoint=abs(bx'*z-x'*btz)/max([1,norm(bx)*norm(z),norm(btz)*norm(x)]);
assert(adjoint<=1e-8,'AugmentedQR:Adjoint','完整增广QR伴随检查失败。');clear x z bx btz;
before=counts;
[y,flag,nativeRelative,iterations,history,normalHistory]=lsqr(@operate,augmentedRhs,1e-12,options.maximumIterations);
nativeCounts=counts;for key=fieldnames(counts)',nativeCounts.(key{1})=counts.(key{1})-before.(key{1});end
coordinates=reshape(y,shape);counts.P=counts.P+1;value=initial+applyRocketAugmentedPatchQR(factor,coordinates);
image=original(value);weighted=weight(value,'notransp');raw=weighted;
if isfield(metric,'peakRowWeights'),raw=raw./metric.peakRowWeights;end
relative=norm(rhs-image,'fro')/denominator;objective=relative^2+(penalty*norm(weighted,'fro')/denominator)^2;
initialObjective=(norm(rhs-initialImage,'fro')/denominator)^2+(penalty*norm(initialMetric,'fro')/denominator)^2;
assert(all(isfinite(value),'all')&&objective<=initialObjective*(1+1e-10),'AugmentedQR:Objective','同一增广目标未下降。');
report=struct('method','complete-parent-C-fixed-primitive-penalty-augmented-patch-QR-LSQR', ...
    'options',options,'penaltyInOriginalResidualUnits',penalty,'relativeCResidual',relative, ...
    'internalCTargetMet',relative<=options.relativeTolerance,'normalizedAugmentedObjective',objective, ...
    'initialNormalizedAugmentedObjective',initialObjective,'primitiveUpdateNorm',norm(weighted,'fro'), ...
    'primitiveUpdateNormIncludesRowWeights',isfield(metric,'peakRowWeights'), ...
    'maximumDensityFraction',max(abs(raw(:,1))),'maximumVelocityOverSound',max(vecnorm(raw(:,2:4),2,2)), ...
    'maximumTemperatureFraction',max(abs(raw(:,5))),'nativeFlag',flag,'nativeIterations',iterations, ...
    'nativeRelativeAugmentedCorrectionResidual',nativeRelative,'nativeResidualHistory',history, ...
    'nativeNormalResidualHistory',normalHistory,'adjointIdentityRelativeDifference',adjoint, ...
    'operatorApplications',counts,'nativeSolverApplications',nativeCounts,'elapsedSeconds',toc(timer), ...
    'originalParentResidualNormUnchanged',true,'originalNSAcceptanceClaimed',false,'physicalSolution',false, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=original(v),counts.C=counts.C+1;out=applyRocketPackedGraph(packed,v);end
    function out=weight(v,flag)
        if strcmp(flag,'notransp'),counts.W=counts.W+1;else,counts.WT=counts.WT+1;end
        out=applyRocketPrimitiveUpdateMetric(metric,v,flag);
    end
    function out=operate(v,flag)
        if strcmp(flag,'notransp')
            counts.P=counts.P+1;mapped=applyRocketAugmentedPatchQR(factor,reshape(v,shape));
            top=original(mapped);bottom=penalty*weight(mapped,'notransp');out=[top(:);bottom(:)];
        elseif strcmp(flag,'transp')
            counts.CT=counts.CT+1;counts.PT=counts.PT+1;
            back=applyRocketPackedGraphTranspose(packed,reshape(v(1:N),shape))+penalty*weight(reshape(v(N+1:end),shape),'transp');
            mapped=applyRocketAugmentedPatchQR(factor,back,'transp');out=mapped(:);
        else,error('AugmentedQR:Flag','未知增广QR作用方向。');end
    end
end
