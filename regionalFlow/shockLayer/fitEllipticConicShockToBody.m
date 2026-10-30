function fit=fitEllipticConicShockToBody(V,F,upstream,gas,seed,coordinates,input)
%FITELLIPTICCONICSHOCKTOBODY 以原面元法线距离匹配内部薄层生成体形，不拟合旧压力。
% 模块一：delta与正定K五个参数共同变化，原V/F与真实上游保持不变。
% 模块二：有限差分Gauss-Newton/阻尼只用于体形残差，不改物理阈值或裁状态。
% 模块三：达不到显式匹配门槛时如实返回未匹配，不把最小二乘最优称完整外流。
if nargin<7,input=struct;end
maximum=8;if isfield(input,'maximumIterations'),maximum=input.maximumIterations;end
tolerance=1e-5;if isfield(input,'bodyMatchToleranceM'),tolerance=input.bodyMatchToleranceM;end
order=8;if isfield(input,'quadratureOrder'),order=input.quadratureOrder;end
validateattributes(maximum,{'double'},{'real','finite','scalar','integer','positive','<=',30});
validateattributes(tolerance,{'double'},{'real','finite','scalar','positive'});
assert(size(coordinates,1)>5&&size(coordinates,2)==2,'ShockLayer:BodyFit','需要多于五个真实匹配列。');
scale=seed.referenceRadiusM;initial=chol(seed.K*scale,'lower');
parameters=[log(seed.standOffM/scale);log(initial(1,1));log(initial(2,2));initial(2,1);seed.profileB*scale];
damping=1e-6;history=struct('iteration',{},'parameters',{},'rmsM',{},'maximumResidualM',{},'damping',{});
evaluations=0;rejectedTrials=0;
current=evaluate(parameters);first=current;stop='迭代预算耗尽，体形仍未匹配';
for iteration=1:maximum
    history(iteration)=struct('iteration',iteration,'parameters',parameters,'rmsM',norm(current.residualM)/sqrt(numel(current.residualM)), ...
        'maximumResidualM',max(abs(current.residualM)),'damping',damping); %#ok<AGROW>
    fprintf('SHOCK_BODY_FIT iteration=%d rmsM=%.9g maxM=%.9g evaluations=%d\n', ...
        iteration,history(iteration).rmsM,history(iteration).maximumResidualM,evaluations);
    if max(abs(current.residualM))<=tolerance,stop='全部采样列通过显式体形距离门槛';break;end
    J=zeros(numel(current.residualM),5);step=1e-3;
    for j=1:5
        shift=zeros(numel(parameters),1);shift(j)=step;plus=evaluate(parameters+shift);minus=evaluate(parameters-shift);
        J(:,j)=(plus.residualM-minus.residualM)/(2*step*scale);
    end
    residual=current.residualM/scale;accepted=false;
    for attempt=1:6
        delta=-(J'*J+damping*eye(5))\(J'*residual);
        for alpha=2.^(0:-1:-6)
            try
                trial=evaluate(parameters+alpha*delta);
                if norm(trial.residualM)<norm(current.residualM)*(1-1e-10)
                    parameters=parameters+alpha*delta;current=trial;accepted=true;break;
                end
            catch ex
                if ~(startsWith(ex.identifier,'ShockLayer:')||startsWith(ex.identifier,'SurfaceFlow:')),rethrow(ex);end
                rejectedTrials=rejectedTrials+1;
            end
        end
        if accepted,damping=max(damping/5,1e-12);break;else,damping=damping*10;end
    end
    if ~accepted,stop='给定五参数激波族没有进一步降低真实体形残差，未放宽门槛';break;end
end
matched=max(abs(current.residualM))<=tolerance;
if matched,stop='全部采样列通过显式体形距离门槛';end
fit=struct('method','共同椭圆截面二次曲线激波五参数内部质量/熵/压力/厚度反算 + 原三角面距离匹配', ...
    'initial',first,'result',current,'history',history,'parameters',parameters, ...
    'bodyMatchToleranceM',tolerance,'shapeMatchingConverged',matched,'stopReason',stop, ...
    'evaluations',evaluations,'rejectedTrials',rejectedTrials,'originalFaceCount',size(F,1), ...
    'scope','五参数薄激波层区域原型；没有驻点极限/跨声速整体及下游推进验收', ...
    'stagnationGlobalMatchingSolved',false,'fullExternalFlowSolved',false, ...
    'boundaryLayerEdgeAvailable',false,'heatFluxSolved',false,'wallTemperatureSolved',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0);

    function value=evaluate(p)
        evaluations=evaluations+1;L=[exp(p(2)) 0;p(4) exp(p(3))];
        shape=struct('frame',seed.frame,'K',(L*L')/scale, ...
            'origin',seed.noseCenter-(scale+scale*exp(p(1)))*seed.frame(:,1),'profileB',p(5)/scale);
        count=size(coordinates,1);columns=cell(count,1);points=zeros(count,3);normals=points;thickness=zeros(count,1);pressure=thickness;
        for k=1:count
            columns{k}=solveEllipticConicShockLayerColumn(upstream,gas,shape,coordinates(k,:)',[1;0], ...
                struct('quadratureOrder',order));
            points(k,:)=columns{k}.shock.position';normals(k,:)=columns{k}.shock.normal';
            thickness(k)=columns{k}.normalDistanceM(end);pressure(k)=columns{k}.state.pressure(end);
        end
        hits=intersectShockNormalsWithBody(V,F,points,normals);
        assert(all(hits.hit),'ShockLayer:BodyFit','候选共同激波没有全部位于原几何外侧。');
        value=struct('shape',shape,'columns',{columns},'geometricMapping',hits, ...
            'predictedNormalDistanceM',thickness,'residualM',thickness-hits.distanceM, ...
            'candidatePressurePa',pressure,'flowValuesAssignedToProductionBody',false);
    end
end
