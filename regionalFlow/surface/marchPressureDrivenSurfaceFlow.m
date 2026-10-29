function result=marchPressureDrivenSurfaceFlow(problem,inputOptions)
%MARCHPRESSUREDRIVENSURFACEFLOW 压力驱动曲面射线族及相邻流线横向度量。
% 模块一：调用者明确提供光滑曲面、同一物性的压力/外缘场及入口方向。
% 模块二：五条相邻射线以固定步长RK4积分曲面切向Euler方程，无Dijkstra选路。
% 模块三：两组标签间距的中心差分投影至n×t得到横向h；不能用总距离
% 替代横向度量。两组差值是分辨率指标，不是严格误差上界，仍需外部收敛核对。
% 模块四：遇物性/接口/积分失效或流管聚焦后保留NaN与失败原因，不补值。
% 范围：有明确入口的光滑附着表面流线；压力与热状态仍为输入，不是
% 完整外流连续方程、激波层、边界层、热流或分离再附求解器。
if nargin<2, inputOptions=struct; end
options=pressureSurfaceFlowOptions(inputOptions);
assert(isstruct(problem)&&isscalar(problem)&&all(isfield(problem, ...
    {'sample','seed','labelCenter','labelStep','labelUnit','arcStationsM','gas','sources'})), ...
    'SurfaceFlow:Problem','缺少采样、入口、标签、弧长、物性或来源。');
assert(isa(problem.sample,'function_handle')&&isa(problem.seed,'function_handle'));
for key={'geometry','pressure','thermodynamics','entryDirection'}
    name=key{1};assert(isfield(problem.sources,name)&& ...
        (ischar(problem.sources.(name))||isstring(problem.sources.(name)))&& ...
        strlength(string(problem.sources.(name)))>0,'SurfaceFlow:Source','每项输入必须明确来源。');
end
assert(any(strcmp(problem.labelUnit,{'m','rad','dimensionless'})),'SurfaceFlow:Label','标签单位须明确。');
validateattributes(problem.labelCenter,{'double'},{'real','finite','scalar'});
validateattributes(problem.labelStep,{'double'},{'real','finite','scalar','positive'});
stations=problem.arcStationsM(:);validateattributes(stations,{'double'},{'real','finite','column','nonnegative'});
assert(numel(stations)>=2&&stations(1)==0&&all(diff(stations)>0),'SurfaceFlow:Stations','弧长须从0严格递增。');
g=problem.gas;
assert(all(isfield(g,{'gamma','R','temperatureRangeK','referenceTemperatureK','referencePressurePa'})));
validateattributes(g.temperatureRangeK,{'double'},{'real','finite','size',[1 2],'positive'});
assert(g.temperatureRangeK(2)>g.temperatureRangeK(1));
validateattributes(g.referenceTemperatureK,{'double'},{'real','finite','scalar','positive'});
validateattributes(g.referencePressurePa,{'double'},{'real','finite','scalar','positive'});
labels=problem.labelCenter+[-1 -0.5 0 0.5 1]*problem.labelStep;rays=cell(5,1);
for ray=1:5
    seed=problem.seed(labels(ray));
    assert(isstruct(seed)&&all(isfield(seed,{'coordinate','direction'})),'SurfaceFlow:Seed','入口须提供坐标与三维方向。');
    validateattributes(seed.coordinate,{'double'},{'real','finite','size',[2 1]});
    validateattributes(seed.direction,{'double'},{'real','finite','size',[3 1]});
    initial=problem.sample(seed.coordinate);A=initial.firstDerivatives;
    lengths=vecnorm(A,2,1)';unitA=A./lengths';unitDirection=seed.direction/norm(seed.direction);
    coordinateDirection=((unitA'*unitA)\(unitA'*unitDirection))./lengths;
    assert(norm(A*coordinateDirection-unitDirection)<=options.maximumUnitTangentError, ...
        'SurfaceFlow:Seed','入口方向必须为曲面切向；不靠投影修正法向输入。');
    state=[seed.coordinate;coordinateDirection];[~,initialInfo]=pressureSurfaceFlowRHS(state,problem,options);
    localProblem=problem;localProblem.referenceTotalEnthalpy=initialInfo.totalEnthalpy;
    localProblem.referenceEntropy=initialInfo.entropy;
    rays{ray}=integrate(state,localProblem,stations,options);rays{ray}.label=labels(ray);
end

%% 相邻标签差分得到横向尺度；聚焦或缺少邻线后不再返回有效h
n=numel(stations);metric=nan(n,1);coarseMetric=metric;metricDifference=metric;
fullSeparation=metric;valid=true(n,1);reason=repmat({''},n,1);
for k=1:n
    if ~all(cellfun(@(r)r.valid(k),rays))
        valid(k:end)=false;reason(k:end)={'相邻射线未完成同一站点；不得外推流管'};break;
    end
    separation=(rays{4}.position(k,:)-rays{2}.position(k,:))/problem.labelStep;
    coarse=(rays{5}.position(k,:)-rays{1}.position(k,:))/(2*problem.labelStep);
    transverse=cross(rays{3}.normal(k,:),rays{3}.tangent(k,:));
    metric(k)=dot(transverse,separation);fullSeparation(k)=norm(separation);
    coarseMetric(k)=dot(transverse,coarse);metricDifference(k)=abs(coarseMetric(k)-metric(k))/3;
    if k==1
        assert(isfinite(metric(k))&&metric(k)>0&&metricDifference(k)<=options.maximumRelativeMetricDifference*metric(k), ...
            'SurfaceFlow:MetricSeed','入口标签须沿n×t正向、宽度非零且两组标签间距相容。');
    elseif ~isfinite(metric(k))||metric(k)<=options.minimumMetricRatio*metric(1)
        valid(k:end)=false;reason(k:end)={'流管聚焦、标签反向或低于声明宽度门槛'};metric(k:end)=NaN;break;
    elseif metricDifference(k)>options.maximumRelativeMetricDifference*metric(k)
        valid(k:end)=false;reason(k:end)={'两组标签间距不能分辨当地横向宽度；不把残余宽度当作有效流管'};metric(k:end)=NaN;break;
    end
end
metric(~valid)=NaN;center=rays{3};
for key={'coordinate','position','tangent','normal','pressure','temperature','density','speed','totalEnthalpy','entropy'}
    name=key{1};center.(name)(~valid,:)=NaN;
end
center.valid=valid;center.missingReason=reason;
result=struct('method','曲面切向Euler压力转向 + 五条RK4射线 + 双间距标签中心差分横向度量', ...
    'scope','声明的光滑附着区、定比热相容外缘输入；不解未知压力、激波层、边界层或分离再附', ...
    'sourceURL','https://ntrs.nasa.gov/api/citations/20040121205/downloads/20040121205.pdf', ...
    'sourceEquations','印刷第4页式13的曲面切向Euler形式；h由相邻射线的横向分离独立构造，非完整VII复现', ...
    'arcLengthM',stations,'center',center,'rays',{rays},'metric',metric, ...
    'metricRatio',metric/metric(1),'coarseMetricDiagnostic',coarseMetric, ...
    'metricDifferenceIndicator',metricDifference,'fullLabelSeparationNorm',fullSeparation, ...
    'centralRayIndexOneBased',3,'fineLabelStep',problem.labelStep/2, ...
    'metricUnit',['m/' problem.labelUnit],'labelStep',problem.labelStep,'labelCenter',problem.labelCenter, ...
    'validStation',valid,'missingReason',{reason},'options',options,'gas',g,'sources',problem.sources, ...
    'pressureStatus','输入场量，未由本模块求解','thermodynamicStatus','输入外缘及状态方程核对，不是新求温度场', ...
    'directionStatus','条件于声明输入场与入口的数值流向','metricStatus','相邻射线标签中心差分，需标签/步长收敛核对', ...
    'rayDiagnosticScope','rays为单独射线诊断；validStation为流管可用前缀，不能用rays补齐失效站点', ...
    'massContinuityOfFullExternalFlowSolved',false,'separationBoundarySolved',false, ...
    'metricDifferenceIsRigorousBound',false,'entryValidationRHSCalls',5, ...
    'physicalStateClipped',false,'productionQualified',false,'newProductionSolvedAreaM2',0);
end

function ray=integrate(state,problem,stations,options)
% 模块用途：每条射线保存实际RK4端点、采样导数和停止原因，拒绝后不重播种。
n=numel(stations);ray=struct('valid',false(n,1),'coordinate',nan(n,2),'position',nan(n,3), ...
    'tangent',nan(n,3),'normal',nan(n,3),'pressure',nan(n,1),'temperature',nan(n,1), ...
    'density',nan(n,1),'speed',nan(n,1),'totalEnthalpy',nan(n,1),'entropy',nan(n,1), ...
    'state',nan(n,4),'unitTangentError',nan(n,1),'samples',{cell(n,1)}, ...
    'stopReason','','failureIdentifier','','failedRequestedStation',NaN,'substeps',0,'RHSCalls',0);
for k=1:n
    try
        if k>1
            count=ceil((stations(k)-stations(k-1))/options.maximumArcStepM);
            assert(ray.substeps+count<=options.maximumSubstepsPerRay,'SurfaceFlow:Budget','超出已声明的射线步数预算。');
            h=(stations(k)-stations(k-1))/count;
            for j=1:count
                a=callRHS(state);b=callRHS(state+h*a/2);
                c=callRHS(state+h*b/2);d=callRHS(state+h*c);
                state=state+h*(a+2*b+2*c+d)/6;ray.substeps=ray.substeps+1;
            end
        end
        [~,info]=callRHS(state);
        assert(info.unitTangentError<=options.maximumUnitTangentError,'SurfaceFlow:ArcLength','弧长切向量漂移超限；不重归一积分状态。');
        ray.valid(k)=true;ray.coordinate(k,:)=state(1:2)';ray.state(k,:)=state';
        for name={'position','tangent','normal','totalEnthalpy','entropy','unitTangentError'}
            key=name{1};ray.(key)(k,:)=info.(key)';
        end
        for name={'pressure','temperature','density','speed'}
            key=name{1};ray.(key)(k)=info.sample.(key);
        end
        ray.samples{k}=info.sample;
    catch ex
        if ~startsWith(ex.identifier,'SurfaceFlow:'), rethrow(ex); end
        if k==1, rethrow(ex); end
        ray.stopReason=ex.message;ray.failureIdentifier=ex.identifier;ray.failedRequestedStation=k;return;
    end
end
ray.stopReason='全部请求弧长站点已积分；不代表完整外流或热流已解';
    function [value,details]=callRHS(y)
        ray.RHSCalls=ray.RHSCalls+1;
        [value,details]=pressureSurfaceFlowRHS(y,problem,options);
    end
end
