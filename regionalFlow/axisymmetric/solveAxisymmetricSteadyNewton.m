function result=solveAxisymmetricSteadyNewton(Q,mesh,gas,boundary,options)
%SOLVEAXISYMMETRICSTEADYNEWTON 直接求原NS/RANS空间残差为零的开放域定常Newton解。
% 模块一：四/五量原始空间方程，不含伪时间质量项，不改变物理边界或温域。
% 模块二：局部尺度完整有色Jacobian、稀疏Newton与原型Armijo回溯；不裁状态。
% 模块三：原2e-6空间门槛及边界/源通量账同时检查；不是BE接受或物理/网格资格。
% 闭域需要额外守恒约束，本切片明确拒绝；不靠任意钉点建立假唯一解。
if nargin<5,options=struct;end
defaults=struct('maximumIterations',24,'maximumLineSearch',14,'residualScale',[], ...
    'steadyResidualTolerance',2e-6,'conservationTolerance',5e-10,'differenceStep',2e-6);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'AxisRoot:Options','未知直接定常选项。');
for k=fieldnames(defaults)',if ~isfield(options,k{1}),options.(k{1})=defaults.(k{1});end,end
for k={'maximumIterations','maximumLineSearch'},validateattributes(options.(k{1}),{'numeric'},{'real','finite','positive','integer','scalar'});end
for k={'steadyResidualTolerance','conservationTolerance','differenceStep'},validateattributes(options.(k{1}),{'numeric'},{'real','finite','positive','scalar'});end
nq=size(Q,2);assert(ismember(nq,[4 5])&&size(Q,1)==mesh.cellCount,'AxisRoot:State','需要匹配的四/五量状态。');
validateattributes(options.residualScale,{'numeric'},{'real','finite','positive','numel',nq});scale=reshape(options.residualScale,1,nq);
open=false;for key={'iMin','iMax','jMin','jMax'},open=open||~ismember(boundary.(key{1}).type,{'wall','axis','periodic'});end
assert(open,'AxisRoot:ClosedDomain','本切片仅支持明确开放域，闭域不得省略守恒约束。');
initialQ=Q;U=Q./scale;[rate,info]=evaluate(Q);evaluations=1;norm=max(abs(rate)./scale,[],'all');balance=steadyBalance(info);
history=struct('iteration',0,'residual',norm,'boundaryBalance',balance,'alpha',0,'linearDefect',0,'residualEvaluations',1);
pattern=[];converged=false;reason='maximum-iterations';accepted=0;linearSolves=0;jacobians=0;lineLog=struct([]);
for iteration=1:options.maximumIterations
    if norm<=options.steadyResidualTolerance&&balance<=options.conservationTolerance,converged=true;reason='original-spatial-root-and-flux-balance-met';break;end
    g=rate./scale;
    [J,jreport,pattern]=axisymmetricFVJacobianLocal(Q,mesh,gas,boundary,scale, ...
        struct('differenceStep',options.differenceStep,'differenceScheme','central'),pattern);
    evaluations=evaluations+jreport.residualEvaluations;jacobians=jacobians+1;
    if ~jreport.success,reason='jacobian-perturbation-failed';break;end
    du=reshape(J\(-g(:)),size(U));linearSolves=linearSolves+1;
    linearDefect=max(abs(J*du(:)+g(:)))/max(1,norm);
    if ~all(isfinite(du),'all'),reason='nonfinite-linear-update';break;end
    found=false;alpha=1;
    for search=0:options.maximumLineSearch
        identifier='';candidateNorm=NaN;candidateBalance=NaN;evaluations=evaluations+1;
        try
            candidate=U+alpha*du;[newRate,newInfo]=evaluate(candidate.*scale);
            candidateNorm=max(abs(newRate)./scale,[],'all');candidateBalance=steadyBalance(newInfo);
            found=candidateNorm<=(1-1e-4*alpha)*norm|| ...
                (candidateNorm<=options.steadyResidualTolerance&&candidateBalance<=options.conservationTolerance);
        catch exception
            identifier=exception.identifier;
            if ~ismember(identifier,{'AxisFV:State','AxisFV:Temperature','AxisSA:State','AxisFV:HLLC','AxisFV:HLLE', ...
                    'AxisFV:PressureSource','AxisFV:AxisTraceState','AxisSA:AxisTraceState','AxisFV:SAModelState', ...
                    'AxisFV:TransportTemperature','AxisFV:TransportState','AxisFV:WallThermalState','AxisSA:Arithmetic'}),rethrow(exception);end
        end
        row=struct('iteration',iteration,'search',search,'alpha',alpha,'baseNorm',norm,'candidateNorm',candidateNorm, ...
            'candidateBoundaryBalance',candidateBalance,'accepted',found,'errorIdentifier',identifier);
        if isempty(lineLog),lineLog=row;else,lineLog(end+1)=row;end %#ok<AGROW>
        if found,break;end;alpha=alpha/2;
    end
    if ~found,reason='line-search-failed';break;end
    U=candidate;Q=U.*scale;rate=newRate;info=newInfo;norm=candidateNorm;balance=candidateBalance;accepted=accepted+1;
    history.iteration(end+1)=accepted;history.residual(end+1)=norm;history.boundaryBalance(end+1)=balance;
    history.alpha(end+1)=alpha;history.linearDefect(end+1)=linearDefect;history.residualEvaluations(end+1)=evaluations;
    fprintf('AxisRoot iteration=%d original residual=%.9g alpha=%.6g flux balance=%.3g calls=%d\n',accepted,norm,alpha,balance,evaluations);
end
if norm<=options.steadyResidualTolerance&&balance<=options.conservationTolerance,converged=true;reason='original-spatial-root-and-flux-balance-met';end
result=struct('kind','open-axisymmetric-original-spatial-Newton-root-v1','initialQ',initialQ,'Q',Q, ...
    'finalRate',rate,'finalInfo',info,'residualNorm',norm,'boundaryFluxBalanceRelative',balance, ...
    'spatialRootConverged',converged,'spatialRootToleranceMet',norm<=options.steadyResidualTolerance, ...
    'acceptedNewtonIterations',accepted,'residualEvaluations',evaluations,'linearSolves',linearSolves, ...
    'jacobianBuilds',jacobians,'history',history,'lineLog',lineLog,'stopReason',reason,'options',options, ...
    'residualScale',scale,'physicalTimeAdvanced',0,'BETimeStepsClaimed',0,'originalSteadyToleranceChanged',false, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
    function [r,i]=evaluate(state)
        if nq==5,[r,i]=axisymmetricRANSResidual(state,mesh,gas,boundary);else,[r,i]=axisymmetricNavierStokesResidual(state,mesh,gas,boundary);end
    end
    function value=steadyBalance(i)
        magnitude=sum(abs(i.faceFluxIntegralPerRadian(mesh.boundaryFaces,:)),1)+sum(abs(i.cellSourceIntegral),1);
        net=abs(i.netBoundaryFlux-i.sourceIntegral);ratios=zeros(size(net));nonzero=magnitude>0;
        ratios(nonzero)=net(nonzero)./magnitude(nonzero);ratios(~nonzero&net>0)=Inf;value=max(ratios);
    end
end
