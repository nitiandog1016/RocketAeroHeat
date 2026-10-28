function report=studyRocketFiniteStepPatch(folder,fixture,mode)
%STUDYROCKETFINITESTEPPATCH 同真实夹具、同预算对照有限步一致性保护。
% 先逐量回放原完整NS核心行，再用相同稠密方向预算做连续非线性推进。
% 小域外侧固定为同一Q24邻域；不是整箭边界替代或物理验证。
if nargin<3,mode='agreement-guard';end
assert(ismember(mode,{'agreement-guard','minimum-state','minimum-state-matched'}));
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);s=load(fixture,'patch');patch=s.patch;
[r0,i0]=rocketFiniteStepPatchResidual(patch.Q0,patch);[r1,i1]=rocketFiniteStepPatchResidual(patch.Q1,patch);
replay=[relative(r0,patch.expectedRate0),relative(r1,patch.expectedRate1), ...
    relative(i0.gradient,patch.expectedGradient0),relative(i1.gradient,patch.expectedGradient1), ...
    relative(i0.faceFlux,patch.expectedFlux0),relative(i1.faceFlux,patch.expectedFlux1)];
assert(max(replay)<1e-12,'Patch:Replay','缩小算例不能复现真实NS行/面/重构。');
% 原Q22到Q24路径上只改变步长，原始变量和无量纲表示均回到同一个Q。
direction=patch.direction;h=1e-6/max(abs(direction),[],'all');
[rh,~]=rocketFiniteStepPatchResidual(patch.Q0+h*direction.*patch.scale,patch);jv=(rh-r0)./patch.scale/h;
fractions=[2.^(-12:0)];sweep=struct([]);
for f=fractions
    alpha=f/512;[ra,ia]=rocketFiniteStepPatchResidual(patch.Q0+alpha*direction.*patch.scale,patch);
    delta=(ra-r0)./patch.scale;rem=delta-alpha*jv;
    row=struct('fraction',f,'alpha',alpha,'rms',sqrt(mean((ra./patch.scale).^2,'all')), ...
        'maximum',max(abs(ra./patch.scale),[],'all'),'relativeNonlinearRemainder',norm(rem,'fro')/norm(alpha*jv,'fro'), ...
        'changedFaces',nnz(ia.sampleRegion~=i0.sampleRegion));
    if isempty(sweep),sweep=row;else,sweep(end+1)=row;end %#ok<AGROW>
end
baseline=advance(false);guarded=advance(true);
report=struct('status','tested-exact-interior-replay-and-bounded-local-multistep','replayRelativeErrors',replay, ...
    'coreCells',numel(patch.core),'haloIncludedCells',size(patch.Q0,1),'sameInitialState','Q24', ...
    'sameBudget',struct('maximumIterations',6,'maximumCandidatesPerIteration',17,'denseJacobianColumns',5*numel(patch.core)), ...
    'baseline',baseline.summary,'guarded',guarded.summary,'comparisonMode',mode,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'patch_comparison.mat'),'report','baseline','guarded','sweep','replay','-v7');
fid=fopen(fullfile(folder,'patch_comparison.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
    function result=advance(useGuard)
        Q=patch.Q1;timer=tic;calls=0;rows=struct([]);trials=struct([]);directions={};stop='iteration-budget';
        [rate,info]=eval(Q);rows=stateRow(0,0,NaN,NaN,NaN);scale=patch.scale;core=patch.core;
        for iteration=1:6
            r=rate./scale;J=zeros(numel(r));
            % 原空间算子逐列有限差分：相同预算，不改预条件器/线性迭代次数。
            for column=1:numel(r)
                [a,b]=ind2sub(size(r),column);step=1e-6*max(1,abs(Q(core(a),b)/scale(b)));
                plus=Q;minus=Q;plus(core(a),b)=plus(core(a),b)+step*scale(b);minus(core(a),b)=minus(core(a),b)-step*scale(b);
                rp=eval(plus);rm=eval(minus);col=(rp-rm)./scale/(2*step);J(:,column)=col(:);
            end
            du=reshape(J\(-r(:)),size(r));
            if startsWith(mode,'minimum-state')
                % 两组使用完全相同的列/Jacobian及SVD预算，区别仅在所选方向。
                metric=buildRocketPrimitiveUpdateMetric(Q(core,:),scale,patch.gas);Winv=zeros(numel(r));
                for column=1:numel(r)
                    e=zeros(size(r));e(column)=1;v=invertRocketPrimitiveUpdateMetric(metric,e);Winv(:,column)=v(:);
                end
                [z,selection]=solveRocketMinimumStateCoefficients(J*Winv,-r(:),.0999);
                if strcmp(mode,'minimum-state-matched')
                    % 相同0.0999约束的保守变量范数对照，区分不精确性与状态度量的影响。
                    [conservativeChoice,conservativeSelection]=solveRocketMinimumStateCoefficients(J,-r(:),.0999);
                    assert(conservativeSelection.feasible);du=reshape(conservativeChoice,size(r));
                end
                if useGuard
                    if ~selection.feasible,stop='minimum-state-span-infeasible';break;end
                    du=reshape(Winv*z,size(r));
                end
            end
            linearDefect=norm(J*du(:)+r(:))/norm(r(:));
            full=zeros(size(Q));full(core,:)=du;step=1e-6/max(abs(full),[],'all');
            rd=eval(Q+step*full.*scale);j=(rd-rate)./scale/step;
            nsDefect=norm(j+r,'fro')/norm(r,'fro');directions{end+1}=struct('Q',Q,'direction',full,'jv',j,'linearDefect',linearDefect,'nsDefect',nsDefect); %#ok<AGROW>
            if nsDefect>.1,stop='original-direction-gate';break;end
            baseRms=sqrt(mean(r.^2,'all'));baseMax=max(abs(r),[],'all');accepted=false;
            for search=0:16
                alpha=2^-search;err='';rr=NaN;mm=NaN;ratio=NaN;
                try
                    candidate=Q+alpha*full.*scale;[next,ni]=eval(candidate);nextU=next./scale;
                    rr=sqrt(mean(nextU.^2,'all'));mm=max(abs(nextU),[],'all');
                    evidence=rocketFiniteStepAgreement(rate,next,j,scale,alpha);
                    ratio=evidence.relativeNonlinearRemainder;
                    accepted=rr<=(1-1e-4*alpha)*baseRms&&ni.spatialConservationRelative<=1e-9;
                    if useGuard&&strcmp(mode,'agreement-guard'),accepted=accepted&&evidence.acceptable;end
                catch exception
                    if ~ismember(exception.identifier,{'CutNS:State','CutNS:InternalEnergy','CutNS:Temperature', ...
                            'CutExactEuler:State','CutExactEuler:InternalEnergy','CutExactEuler:Temperature','CutExactEuler:Vacuum'}),rethrow(exception);end
                    err=exception.identifier;
                end
                row=struct('iteration',iteration,'alpha',alpha,'rms',rr,'maximum',mm,'relativeNonlinearRemainder',ratio,'accepted',accepted,'errorIdentifier',err);
                if isempty(trials),trials=row;else,trials(end+1)=row;end %#ok<AGROW>
                if accepted,break;end
            end
            if ~accepted,stop='candidate-budget';break;end
            Q=candidate;rate=next;info=ni;rows(end+1)=stateRow(iteration,alpha,linearDefect,nsDefect,ratio); %#ok<AGROW>
            fprintf('PATCH guard=%d iteration=%d rms=%.9g max=%.9g alpha=%.9g remainder=%.6g\n',useGuard,iteration,rr,mm,alpha,ratio);
            if mm<=2e-6&&info.patchSteadyFluxRatio<=5e-10,stop='local-original-tolerances';break;end
        end
        summary=struct('history',rows,'acceptedSteps',numel(rows)-1,'stopReason',stop,'residualEvaluations',calls, ...
            'elapsedSeconds',toc(timer),'finiteStepGuard',useGuard&&strcmp(mode,'agreement-guard'), ...
            'minimumStateSelection',useGuard&&startsWith(mode,'minimum-state'), ...
            'internalLinearTarget',1e-3,'steadyConverged',rows(end).maximum<=2e-6&&rows(end).boundaryFluxRatio<=5e-10);
        result=struct('summary',summary,'Q',Q,'rate',rate,'info',info,'trials',trials,'directions',{directions});
        function [rout,iout]=eval(qin),calls=calls+1;[rout,iout]=rocketFiniteStepPatchResidual(qin,patch);end
        function row=stateRow(it,alpha,lin,ns,rem)
            row=struct('iteration',it,'alpha',alpha,'rms',sqrt(mean((rate./patch.scale).^2,'all')), ...
                'maximum',max(abs(rate./patch.scale),[],'all'),'linearDefect',lin,'originalNSDirectionDefect',ns, ...
                'relativeNonlinearRemainder',rem,'minimumDensity',info.minimumDensity,'temperatureRange',info.temperatureRange, ...
                'spatialConservationRelative',info.spatialConservationRelative,'boundaryFluxRatio',info.patchSteadyFluxRatio,'elapsedSeconds',toc(timer));
        end
    end
end
function value=relative(a,b),value=norm(a(:)-b(:))/max(norm(b(:)),realmin);end
