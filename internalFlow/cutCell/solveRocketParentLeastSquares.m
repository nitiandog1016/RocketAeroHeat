function [value,report,coordinates]=solveRocketParentLeastSquares(packed,factor,rhs,initial,options,columnScale)
%SOLVEROCKETPARENTLEASTSQUARES 直接求完整父C的有界全局LSQR，原CV五量/范数不变。
% 模块一：B=C*M*D，伴随=D*M转置*C转置；D可空，完全绕开旧Schur与固定C代表。
% 模块二：保留给定初猜，只求校正；原C的1e-3以内目标按原完整右端验收。
% 模块三：原生flag或法方程提前停止不代表完整C合格，更不能代表原NS/区域物理解。
% 模块四：第三输出是未乘D的原生坐标，只用于MAT证据，不塞入JSON摘要。
% 方法来源：Paige和Saunders的LSQR，ACM TOMS 8(1),1982，doi:10.1145/355984.355989。
if nargin<5,options=struct;end
if nargin<6,columnScale=[];end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'adjointIdentityTolerance',1e-8,'maximumVectorEstimateBytes',256*1024^2);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'ParentLSQR:Options','未知直接父矩阵最小二乘选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});
assert(isscalar(options.relativeTolerance)&&isfinite(options.relativeTolerance)&&options.relativeTolerance>0&&options.relativeTolerance<=1e-3,'ParentLSQR:Options','完整父C目标不得放宽超过1e-3。');
assert(isscalar(options.adjointIdentityTolerance)&&isfinite(options.adjointIdentityTolerance)&&options.adjointIdentityTolerance>0&&options.adjointIdentityTolerance<=1e-8,'ParentLSQR:Options','伴随检查不得放宽超过1e-8。');
n=packed.cellCount;assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&&factor.cellCount==n&& ...
    strcmp(factor.kind,'fixed-five-conservation-strong-connection-local-block-inverse')&&isequal(size(rhs),[n 5])&&isequal(size(initial),size(rhs))&& ...
    isreal(rhs)&&isreal(initial)&&all(isfinite(rhs),'all')&&all(isfinite(initial),'all'),'ParentLSQR:Input','完整父矩阵、原局部映射、右端或初猜非法。');
assert(isempty(columnScale)||(isequal(size(columnScale),size(rhs))&&isreal(columnScale)&&all(isfinite(columnScale)&columnScale>0,'all')),'ParentLSQR:Scale','右侧尺度必须为全部原CV五量的严格正有限数组。');
% 这是声明的20个全向量额外数组估计，不是原生库或进程峰值内存保证。
vectorEstimate=8*numel(rhs)*20;assert(vectorEstimate<=options.maximumVectorEstimateBytes,'ParentLSQR:Budget','声明的额外全向量数组估计超预算。');
counts=struct('C',0,'CT',0,'M',0,'MT',0);solverCounts=counts;timer=tic;
image=original(initial);remaining=rhs-image;denominator=max(realmin,norm(rhs,'fro'));target=options.relativeTolerance*denominator;
initialRelative=norm(remaining,'fro')/denominator;coordinates=zeros(size(rhs));adjointDifference=NaN;reconstruction=0;
nativeFlag=0;nativeRelative=0;iterations=0;nativeTolerance=NaN;residualHistory=norm(remaining,'fro');leastSquaresHistory=[];
if norm(remaining,'fro')>target
    v=sin((1:numel(rhs))');w=cos((1:numel(rhs))'/sqrt(2));bv=operate(v,'notransp');btw=operate(w,'transp');
    adjointDifference=abs(bv'*w-v'*btw)/max([1,norm(bv)*norm(w),norm(btw)*norm(v)]);
    assert(adjointDifference<=options.adjointIdentityTolerance,'ParentLSQR:Adjoint','实际完整父矩阵右预条件伴随检查失败。');clear v w bv btw;
    before=counts;nativeTolerance=min(1e-12,target/norm(remaining,'fro'));
    [y,nativeFlag,nativeRelative,iterations,residualHistory,leastSquaresHistory]=lsqr(@operate,remaining(:),nativeTolerance,options.maximumIterations);
    for key=fieldnames(counts)',solverCounts.(key{1})=counts.(key{1})-before.(key{1});end
    coordinates=reshape(y,size(rhs));q=coordinates;if ~isempty(columnScale),q=q.*columnScale;end
    counts.M=counts.M+1;correction=solveRocketCoupledPatches(factor,q);correctionImage=original(correction);value=initial+correction;image=original(value);
    reconstruction=norm((rhs-image)-(remaining-correctionImage),'fro')/max(1,denominator);
    assert(all(isfinite(value),'all')&&reconstruction<1e-8,'ParentLSQR:Reconstruction','完整父C总方向与校正像不一致。');
else,value=initial;end
relative=norm(rhs-image,'fro')/denominator;
report=struct('method','direct-complete-parent-C-fixed-local-map-right-image-scale-LSQR','options',options, ...
    'initialRelativeCResidual',initialRelative,'relativeCResidual',relative,'internalCTargetMet',relative<=options.relativeTolerance, ...
    'nativeFlag',nativeFlag,'nativeRelativeCorrectionResidual',nativeRelative,'nativeIterations',iterations,'nativeToleranceRelativeToCorrectionRHS',nativeTolerance, ...
    'nativeResidualHistory',residualHistory,'nativeLeastSquaresHistory',leastSquaresHistory,'adjointIdentityRelativeDifference',adjointDifference, ...
    'totalDirectionReconstructionDifference',reconstruction,'operatorApplications',counts,'nativeSolverApplications',solverCounts, ...
    'columnScalingUsed',~isempty(columnScale),'declaredExtraVectorArrayEstimateBytes',vectorEstimate,'processPeakMemoryClaimed',false,'elapsedSeconds',toc(timer), ...
    'completeParentOperatorFixed',true,'baseLocalMapFixedLinear',true,'inverseFromBoundedLSQRClaimedFixedLinear',false, ...
    'oldSchurUsed',false,'fixedCRepresentativesUsed',false,'originalParentResidualNormUnchanged',true, ...
    'independentAllRealBidiagonalizationColumnsClaimed',false,'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=original(update),counts.C=counts.C+1;out=applyRocketPackedGraph(packed,update);end
    function out=operate(v,flag)
        if strcmp(flag,'notransp')
            q=reshape(v,size(rhs));if ~isempty(columnScale),q=q.*columnScale;end
            counts.M=counts.M+1;mapped=solveRocketCoupledPatches(factor,q);result=original(mapped);
        elseif strcmp(flag,'transp')
            counts.CT=counts.CT+1;counts.MT=counts.MT+1;mapped=applyRocketPackedGraphTranspose(packed,reshape(v,size(rhs)));
            result=solveRocketCoupledPatchesTranspose(factor,mapped);if ~isempty(columnScale),result=result.*columnScale;end
        else,error('ParentLSQR:TransposeFlag','未知完整父矩阵LSQR作用方向。');end
        out=result(:);
    end
end
