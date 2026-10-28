function [value,report,coordinates]=solveRocketTemperatureProjectedDirection(packed,factor,rhs,metric,projector,rows,maximumIterations)
%SOLVEROCKETTEMPERATUREPROJECTEDDIRECTION 在完整增广目标上求采样温度一阶不变的研究方向。
% 模块一：v=P*Pi*y，优化原C与冻结W/penalty，显式增加H*v=0方向限制。
% 模块二：正向C*P*Pi，转置Pi*P转置*C转置；不把投影放在错误一侧。
% 模块三：原C1e-3、完整模型/NS0.1及候选温域分别验收，投影不是物理求解或温度裁剪。
if nargin<7,maximumIterations=128;end
validateattributes(maximumIterations,{'numeric'},{'scalar','integer','positive'});
shape=[packed.cellCount 5];N=prod(shape);penalty=factor.penaltyInOriginalResidualUnits;
assert(isequal(size(rhs),shape)&&metric.cellCount==packed.cellCount&&projector.cellCount==packed.cellCount&& ...
    rows.cellCount==packed.cellCount&&all(isfinite(rhs),'all')&&8*N*32<=384*1024^2, ...
    'TemperatureProjection:Input','当前增广系统或向量预算非法。');
timer=tic;counts=struct('C',0,'CT',0,'P',0,'PT',0,'Pi',0,'W',0,'WT',0);
denominator=max(realmin,norm(rhs,'fro'));augmentedRhs=[rhs(:);zeros(N,1)];
x=sin((1:N)');z=cos((1:2*N)'/sqrt(2));bx=operate(x,'notransp');btz=operate(z,'transp');
adjoint=abs(bx'*z-x'*btz)/max([1,norm(bx)*norm(z),norm(btz)*norm(x)]);
assert(adjoint<=1e-8,'TemperatureProjection:Adjoint','增广温度投影完整伴随检查失败。');clear x z bx btz;
before=counts;
[y,flag,nativeRelative,iterations,history,normalHistory]=lsqr(@operate,augmentedRhs,1e-12,maximumIterations);
nativeCounts=counts;for key=fieldnames(counts)',nativeCounts.(key{1})=counts.(key{1})-before.(key{1});end
coordinates=reshape(y,shape);projected=applyRocketTemperatureProjector(projector,coordinates);
value=applyRocketAugmentedPatchQR(factor,projected);image=applyRocketPackedGraph(packed,value);
weighted=applyRocketPrimitiveUpdateMetric(metric,value);raw=weighted;
if isfield(metric,'peakRowWeights'),raw=raw./metric.peakRowWeights;end
relative=norm(rhs-image,'fro')/denominator;objective=relative^2+(penalty*norm(weighted,'fro')/denominator)^2;
derivative=rows.H*value(:);constraintError=norm(derivative)/max(1,norm(rows.H,'fro')*norm(value,'fro'));
assert(all(isfinite(value),'all')&&objective<=1+1e-10&&constraintError<=1e-10, ...
    'TemperatureProjection:Objective','受限方向目标未下降或温度一阶约束未满足。');
report=struct('method','complete-parent-C-fixed-primitive-penalty-temperature-row-nullspace-QR-LSQR', ...
    'maximumIterations',maximumIterations,'penaltyInOriginalResidualUnits',penalty, ...
    'relativeCResidual',relative,'internalCTarget',1e-3,'internalCTargetMet',relative<=1e-3, ...
    'normalizedAugmentedObjective',objective,'primitiveUpdateNorm',norm(weighted,'fro'), ...
    'maximumDensityFraction',max(abs(raw(:,1))),'maximumVelocityOverSound',max(vecnorm(raw(:,2:4),2,2)), ...
    'maximumTemperatureFraction',max(abs(raw(:,5))),'constraintDerivativeK',derivative, ...
    'constraintNormalizedError',constraintError,'constraintFaceIDs',rows.faceIDs, ...
    'nativeFlag',flag,'nativeIterations',iterations,'nativeRelativeAugmentedResidual',nativeRelative, ...
    'nativeResidualHistory',history,'nativeNormalResidualHistory',normalHistory, ...
    'adjointIdentityRelativeDifference',adjoint,'nativeSolverApplications',nativeCounts,'operatorApplications',counts, ...
    'elapsedSeconds',toc(timer),'zeroInitialGuess',true,'physicalBoundaryConditionAdded',false, ...
    'nonlinearTemperatureBoundGuaranteed',false,'constrainedOptimumClaimed',false,'originalNSAcceptanceClaimed',false, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'physicalSolution',false);
    function out=operate(v,flag)
        if strcmp(flag,'notransp')
            counts.Pi=counts.Pi+1;counts.P=counts.P+1;
            projected=applyRocketTemperatureProjector(projector,reshape(v,shape));mapped=applyRocketAugmentedPatchQR(factor,projected);
            counts.C=counts.C+1;counts.W=counts.W+1;top=applyRocketPackedGraph(packed,mapped);bottom=penalty*applyRocketPrimitiveUpdateMetric(metric,mapped);
            out=[top(:);bottom(:)];
        elseif strcmp(flag,'transp')
            counts.CT=counts.CT+1;counts.WT=counts.WT+1;counts.PT=counts.PT+1;counts.Pi=counts.Pi+1;
            back=applyRocketPackedGraphTranspose(packed,reshape(v(1:N),shape))+ ...
                penalty*applyRocketPrimitiveUpdateMetric(metric,reshape(v(N+1:end),shape),'transp');
            mapped=applyRocketAugmentedPatchQR(factor,back,'transp');mapped=applyRocketTemperatureProjector(projector,mapped);out=mapped(:);
        else,error('TemperatureProjection:Flag','未知完整增广投影作用方向。');end
    end
end
