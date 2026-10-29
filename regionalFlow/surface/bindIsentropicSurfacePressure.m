function problem=bindIsentropicSurfacePressure(problem,totalState)
%BINDISENTROPICSURFACEPRESSURE 用固定总状态把连续压力采样接成相容外缘。
% 模块一：调用者提供光滑附着曲面、连续p及其两向导数，来源不能省略。
% 模块二：同一补片只使用一个明确p0/T0；由熵与总焓反算T/rho/U。
% 模块三：由同一状态解析求导，不能逐站重置总压或跨激波/分离接口。
% 范围：单一总状态的等熵表面补片；本函数不求压力，不证明完整外流质量守恒。
assert(isstruct(problem)&&isscalar(problem)&&isa(problem.sample,'function_handle')&& ...
    isfield(problem,'gas')&&isfield(problem,'sources'),'SurfaceFlow:PressureAdapter','缺少压力采样、物性或来源。');
required={'totalPressure','totalTemperature','thermalPropertyMode','method'};
assert(isstruct(totalState)&&isscalar(totalState)&&all(isfield(totalState,required)), ...
    'SurfaceFlow:PressureAdapter','总状态必须有明确温压、物性模式和来源方法。');
if isfield(totalState,'available')
    assert(isequal(totalState.available,true),'SurfaceFlow:PressureAdapter','不能接入求解失败的总状态。');
end
validateattributes(totalState.totalPressure,{'double'},{'real','finite','scalar','positive'});
validateattributes(totalState.totalTemperature,{'double'},{'real','finite','scalar','positive'});
assert(strlength(string(totalState.method))>0,'SurfaceFlow:Source','总状态必须明确来源。');
g=problem.gas;reference=surfaceFlowThermodynamicState(totalState.totalTemperature,totalState.totalPressure,g);
assert(strcmp(reference.thermalPropertyMode,totalState.thermalPropertyMode), ...
    'SurfaceFlow:Gas','压力锚点和流线物性必须一致。');
if strcmp(reference.thermalPropertyMode,'thermallyPerfect')
    assert(isfield(totalState,'thermalPropertyConvention')&& ...
        strcmp(totalState.thermalPropertyConvention,g.thermalPropertyConvention), ...
        'SurfaceFlow:ThermalConvention','旧NASA7总状态必须先从物理锚点显式转换为连续积分约定。');
end
pressureSample=problem.sample;p0=totalState.totalPressure;T0=totalState.totalTemperature;
cachedPressure=NaN;cachedThermal=[];
problem.sample=@sample;
problem.sources.thermodynamics=[char(totalState.method) '；单一p0/T0、同物性熵/总焓恢复与解析导数'];
problem.thermalPressureCoupling=struct('method','单一总状态的等熵压力恢复', ...
    'anchorMethod',totalState.method,'thermalPropertyMode',reference.thermalPropertyMode, ...
    'totalPressure',p0,'totalTemperature',T0,'totalEnthalpy',reference.enthalpy, ...
    'referenceEntropy',reference.entropy,'pressureSolvedByAdapter',false, ...
    'fullExternalContinuitySolved',false,'totalPressureResetAlongRay',false);
if isfield(totalState,'thermalPropertyConvention')
    problem.thermalPressureCoupling.thermalPropertyConvention=totalState.thermalPropertyConvention;
end
if isfield(totalState,'sourceAnchor'),problem.thermalPressureCoupling.anchorConversion=totalState;end

    function value=sample(q)
        value=pressureSample(q);
        assert(isstruct(value)&&isscalar(value)&&all(isfield(value, ...
            {'position','firstDerivatives','secondDerivatives','pressure','pressureGradient','isSmooth','isAttached'})), ...
            'SurfaceFlow:PressureAdapter','缺少压力/曲面及两向导数。');
        assert(isequal(value.isSmooth,true)&&isequal(value.isAttached,true), ...
            'SurfaceFlow:Unsupported','压力补片越过光滑附着区边界。');
        assert(~any(isfield(value,{'temperature','density','speed','totalPressure','totalPressureGradient'})), ...
            'SurfaceFlow:PressureAdapter','压力恢复不能静默覆盖已有外缘或逐站总压；须明确选择输入。');
        p=value.pressure;dp=value.pressureGradient;
        validateattributes(p,{'double'},{'real','finite','scalar','positive'});
        validateattributes(dp,{'double'},{'real','finite','size',[2 1]});
        assert(p<p0,'SurfaceFlow:ThermalPressure','需要0<p<p0的非驻点补片；不裁压力或猜速度。');
        if strcmp(reference.thermalPropertyMode,'thermallyPerfect')
            % 同一闭包中的物性与总状态已冻结；只缓存完全相同的静压，避免直锥重复求根。
            if isequal(p,cachedPressure)
                recovered=cachedThermal;
            else
                recovered=recoverContinuous(p);
                cachedPressure=p;cachedThermal=recovered;
            end
            assert(recovered.valid,'SurfaceFlow:ThermalPressure','同物性压力恢复无解，不回退。');
            T=recovered.temperature;rho=recovered.density;U=recovered.speed;
        else
            T=T0*(p/p0)^((g.gamma-1)/g.gamma);rho=p/(g.R*T);
            U2=2*reference.specificHeat*(T0-T);
            assert(U2>0,'SurfaceFlow:ThermalPressure','外缘动能必须为正。');U=sqrt(U2);
        end
        thermal=surfaceFlowThermodynamicState(T,p,g);cp=thermal.specificHeat;
        assert(U>0,'SurfaceFlow:Stagnation','驻点须独立启动。');
        value.temperature=T;value.density=rho;value.speed=U;
        value.temperatureGradient=g.R*T/(cp*p)*dp;
        value.densityGradient=rho*(dp/p-value.temperatureGradient/T);
        value.speedGradient=-dp/(rho*U);
    end

    function recovered=recoverContinuous(p)
        % 模块用途：对连续标准熵直接二分；与导数中的同一cp(T)严格配套。
        [~,h0,s0]=surfaceFlowAirThermalProperties(T0);
        [~,~,sMin]=surfaceFlowAirThermalProperties(g.temperatureRangeK(1));
        wanted=s0+g.R*log(p/p0);
        assert(wanted>=sMin,'SurfaceFlow:ThermalPressure','压力恢复超出物性下界。');
        lo=g.temperatureRangeK(1);hi=T0;
        for iteration=1:47
            mid=(lo+hi)/2;[~,~,sMid]=surfaceFlowAirThermalProperties(mid);
            if sMid>wanted,hi=mid;else,lo=mid;end
        end
        T=(lo+hi)/2;[~,h,~]=surfaceFlowAirThermalProperties(T);U2=2*(h0-h);
        assert(isfinite(U2)&&U2>0,'SurfaceFlow:ThermalPressure','恢复动能无效，不裁能量。');
        recovered=struct('valid',true,'temperature',T,'density',p/(g.R*T),'speed',sqrt(U2));
    end
end
