function result=solveRocketShellSteadyLine(Q,cut,gas,boundary,options)
%SOLVEROCKETSHELLSTEADYLINE 以五量法向线近似预条件对原完整三维稳态残差做缺陷修正。
% 模块一：B=I/dtau-Japprox仅生成方向；候选必须重新调用原HLLC/WLS/黏性/导热残差。
% 模块二：显式RMS残差下降规则、正性/温域与空间五量通量账；失败候选不提交。
% 模块三：pseudoDt不是物理时间，不称BE方程已解；原最大残差与定常边界账才决定稳态。
if nargin<5,options=struct;end
defaults=struct('maximumIterations',6,'initialPseudoDt',1e-8,'maximumPseudoDt',1e-4, ...
    'minimumPseudoDt',1e-13,'pseudoDtGrowth',3,'maximumStepRetries',3,'maximumLineSearch',8, ...
    'residualScale',[],'steadyResidualTolerance',2e-6,'steadyFluxTolerance',5e-10, ...
    'spatialConservationTolerance',1e-9,'sufficientDecrease',1e-4,'saveTrace',true,'saveFirstOperator',true);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'ShellSteady:Options','未知线预条件迭代选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
for key={'maximumIterations','maximumStepRetries','maximumLineSearch'}
    validateattributes(options.(key{1}),{'numeric'},{'real','finite','integer','scalar','nonnegative'});
end
assert(options.maximumIterations>0);
for key={'initialPseudoDt','maximumPseudoDt','minimumPseudoDt','pseudoDtGrowth','steadyResidualTolerance', ...
        'steadyFluxTolerance','spatialConservationTolerance','sufficientDecrease'}
    validateattributes(options.(key{1}),{'numeric'},{'real','finite','scalar','positive'});
end
assert(options.minimumPseudoDt<=options.initialPseudoDt&&options.initialPseudoDt<=options.maximumPseudoDt&& ...
    options.pseudoDtGrowth>=1&&options.sufficientDecrease<1,'ShellSteady:Options','伪步范围或下降系数不符。');
validateattributes(options.residualScale,{'numeric'},{'real','finite','positive','numel',5});scale=reshape(options.residualScale,1,5);
checkState(Q);initialQ=Q;evaluations=0;[rate,info,spatialBalance,fluxBalance]=evaluate(Q);
initialRate=rate;initialInfo=info;[maximum,merit]=norms(rate);initialMaximum=maximum;initialMerit=merit;
tau=options.initialPseudoDt;accepted=0;trials=0;builds=0;linearSolves=0;trace={};firstOperator=[];
history=struct('iteration',0,'pseudoDt',0,'alpha',0,'maximumResidualPerSecond',maximum,'rmsMerit',merit, ...
    'steadyBoundaryFluxBalance',fluxBalance,'spatialConservationRelative',spatialBalance,'linearDefect',0, ...
    'operatorBuildSeconds',0,'linearSolveSeconds',0,'residualEvaluations',evaluations);
lineLog=struct([]);reason='maximum-nonlinear-iterations';steady=false;
for iteration=1:options.maximumIterations
    if maximum<=options.steadyResidualTolerance&&fluxBalance<=options.steadyFluxTolerance,steady=true;reason='original-spatial-residual-and-steady-boundary-balance-met';break;end
    before=Q;beforeRate=rate;beforeMerit=merit;found=false;
    for retry=0:options.maximumStepRetries
        assert(tau>=options.minimumPseudoDt,'ShellSteady:MinimumPseudoDt','伪步已低于显式下限，没有提交失败状态。');
        ticBuild=tic;operator=buildRocketShellLineOperator(Q,cut,gas,boundary,scale,tau);buildSeconds=toc(ticBuild);builds=builds+1;
        ticSolve=tic;[direction,linear]=solveRocketShellLineSystem(operator,rate./scale);solveSeconds=toc(ticSolve);linearSolves=linearSolves+1;
        if options.saveFirstOperator&&isempty(firstOperator),firstOperator=operator;end
        fprintf('SHELL_LINE_BUILD iteration=%d retry=%d dtau=%.9g build=%.3g solve=%.3g linear=%.3g\n',iteration,retry,tau,buildSeconds,solveSeconds,linear.relativeLinearDefect);
        alpha=1;
        for search=0:options.maximumLineSearch
            trials=trials+1;identifier='';candidateMerit=NaN;candidateMaximum=NaN;candidateFluxBalance=NaN;
            try
                candidate=Q+alpha*direction.*scale;checkState(candidate);
                [newRate,newInfo,newSpatial,newFlux]=evaluate(candidate);[candidateMaximum,candidateMerit]=norms(newRate);candidateFluxBalance=newFlux;
                found=candidateMerit<=(1-options.sufficientDecrease*alpha)*merit;
            catch exception
                identifier=exception.identifier;
                recoverable={'CutNS:State','CutNS:InternalEnergy','CutNS:Temperature','CutNS:HLLC','CutNS:HLLCStarPressure', ...
                    'CutNS:BoundaryMach','RocketAeroHeat:NSDraft:CharacteristicState','RocketAeroHeat:NSDraft:CharacteristicSound', ...
                    'RocketAeroHeat:NSDraft:CharacteristicRegime','RocketAeroHeat:NSDraft:CharacteristicSonic', ...
                    'RocketAeroHeat:NSDraft:CharacteristicAmbiguous'};
                if ~ismember(identifier,recoverable),rethrow(exception);end
            end
            row=struct('iteration',iteration,'retry',retry,'search',search,'pseudoDt',tau,'alpha',alpha, ...
                'baseRmsMerit',merit,'candidateRmsMerit',candidateMerit,'candidateMaximumResidual',candidateMaximum, ...
                'candidateSteadyBoundaryBalance',candidateFluxBalance,'accepted',found,'errorIdentifier',identifier);
            if isempty(lineLog),lineLog=row;else,lineLog(end+1)=row;end %#ok<AGROW>
            fprintf('SHELL_LINE_TRIAL iteration=%d search=%d alpha=%.6g rms=%.9g max=%.9g accepted=%d error=%s\n',iteration,search,alpha,candidateMerit,candidateMaximum,found,identifier);
            if found,break;end;alpha=alpha/2;
        end
        if found,break;end
        tau=tau/4;
    end
    if ~found,reason='bounded-line-search-and-pseudostep-retries-exhausted';break;end
    Q=candidate;rate=newRate;info=newInfo;spatialBalance=newSpatial;fluxBalance=newFlux;maximum=candidateMaximum;merit=candidateMerit;accepted=accepted+1;
    if options.saveTrace
        trace{end+1}=struct('iteration',accepted,'initialQ',before,'initialRate',beforeRate,'linearUpdateU',direction, ...
            'pseudoDt',tau,'alpha',alpha,'Q',Q,'rate',rate,'initialRmsMerit',beforeMerit,'rmsMerit',merit, ...
            'spatialConservationRelative',spatialBalance,'linearReport',linear,'physicalTimeAdvanced',0); %#ok<AGROW>
    end
    for key={'iteration','pseudoDt','alpha','maximumResidualPerSecond','rmsMerit','steadyBoundaryFluxBalance', ...
             'spatialConservationRelative','linearDefect','operatorBuildSeconds','linearSolveSeconds','residualEvaluations'}
        name=key{1};switch name
            case 'iteration',value=accepted;case 'pseudoDt',value=tau;case 'alpha',value=alpha;
            case 'maximumResidualPerSecond',value=maximum;case 'rmsMerit',value=merit;
            case 'steadyBoundaryFluxBalance',value=fluxBalance;case 'spatialConservationRelative',value=spatialBalance;
            case 'linearDefect',value=linear.relativeLinearDefect;case 'operatorBuildSeconds',value=buildSeconds;
            case 'linearSolveSeconds',value=solveSeconds;case 'residualEvaluations',value=evaluations;
        end
        history.(name)(end+1)=value;
    end
    if alpha==1,tau=min(options.maximumPseudoDt,tau*options.pseudoDtGrowth);end
end
if maximum<=options.steadyResidualTolerance&&fluxBalance<=options.steadyFluxTolerance,steady=true;reason='original-spatial-residual-and-steady-boundary-balance-met';end
result=struct('kind','full-three-dimensional-original-NS-line-preconditioned-defect-correction', ...
    'initialQ',initialQ,'initialRate',initialRate,'initialInfo',initialInfo,'Q',Q,'rate',rate,'info',info, ...
    'initialMaximumResidual',initialMaximum,'maximumResidual',maximum,'initialRmsMerit',initialMerit,'rmsMerit',merit, ...
    'acceptedNonlinearIterations',accepted,'candidateTrials',trials,'residualEvaluations',evaluations,'operatorBuilds',builds, ...
    'linearSolves',linearSolves,'history',history,'lineLog',lineLog,'trace',{trace},'firstOperator',firstOperator, ...
    'steadyConverged',steady,'stopReason',reason,'options',options,'residualScale',scale, ...
    'spatialConservationRelative',spatialBalance,'steadyBoundaryFluxBalance',fluxBalance, ...
    'physicalTimeAdvanced',0,'BETimeStepsClaimed',0,'originalSpatialEquationReplaced',false, ...
    'positivityClippingCount',0,'physicalSolution',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
    function checkState(q)
        assert(isequal(size(q),[numel(cut.cells.volume) 5])&&all(isfinite(q),'all')&&isreal(q)&&all(q(:,1)>0),'CutNS:State','候选五量/密度非法。');
        rho=q(:,1);u=q(:,2:4)./rho;internal=q(:,5)-.5*rho.*sum(u.^2,2);assert(all(internal>0),'CutNS:InternalEnergy','候选内能非正，不裁剪。');
        T=internal./(rho*gas.R/(gas.gamma-1));assert(all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)), ...
            'CutNS:Temperature','候选温度超出原物性域，不裁温。');
    end
    function [r,i,spatial,flux]=evaluate(q)
        evaluations=evaluations+1;[r,i]=cutCellNavierStokesResidual(q,cut,gas,boundary);
        f=cut.openFaces;mask=f.cellPlus==0;w=cut.walls;
        magnitude=sum(abs(i.faceFlux(mask,:)).*f.area(mask),1)+sum(abs(i.wallFlux).*w.area,1);
        denom=max(1,magnitude);spatial=max(abs(i.balanceResidual)./denom);flux=max(abs(i.netBoundaryFlux)./denom);
        assert(spatial<=options.spatialConservationTolerance,'ShellSteady:Conservation','原空间五量账未通过，不提交候选。');
    end
    function [maximum,rms]=norms(r)
        z=r./scale;maximum=max(abs(z),[],'all');rms=sqrt(mean(z.^2,'all'));
        assert(isfinite(maximum)&&isfinite(rms),'ShellSteady:Nonfinite','原空间残差非有限。');
    end
end
