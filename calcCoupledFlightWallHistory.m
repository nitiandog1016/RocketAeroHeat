function mesh = calcCoupledFlightWallHistory(mesh,history,material,cfg)
%CALCCOUPLEDFLIGHTWALLHISTORY 变飞行条件下的沿程热边界与材料导热耦合。
% 模块一：沿输入历程线性插值工况，按需建立当前外流/准流线快照。
% 模块二：同一时刻反复用非均匀壁温重积分边界层；材料完整厚度
% 剖面跨时刻接续。自动二分的子步也在自己的时间重算飞行条件。
% 模块三：保存实际工况、算法、失效原因和残差；独立记录固定位置
% 的驻点单点壁温。驻点随攻角移动时不把不同材料点拼成同一热历程。
% 本模型是准定常工程耦合，不解非定常外流、横向导热或分离再附。

%% 1. 核对几何顺序，采用预计算时的气动热配置以保证可复现
input=history.couplingInput;
if ~isequal(mesh.V,input.vertices) || ~isequal(mesh.F,input.faces)
    error('RocketAeroHeat:FlightCoupling:Mesh', ...
        '网格或面元顺序已改变，请重新计算飞行历程。');
end
stationCfg=input.configuration;
initial=cfg.flight.Twall;
if stationCfg.flight.Twall~=initial
    error('RocketAeroHeat:FlightCoupling:Initial','初始壁温已改变，请重新计算飞行历程。');
end
flight=history.flight;
nf=size(mesh.F,1);
cachedCondition=[]; snapshot=[]; currentCondition=[]; flowSnapshot=[];
stationSolves=0; boundaryEvaluations=0;
initialBoundary=evaluate(initial*ones(nf,1),history.time(1));

%% 2. 时间相关回调使自适应子步与外流时刻保持一致
wall=solveCoupledWallConduction(@evaluate,initialBoundary, ...
    history.time,material,initial);
wall.method=['飞行条件分段线性插值 + ' wall.method];
wall.scope=[wall.scope '；外流在同一时刻的壁温迭代内固定，不作粘性-无粘耦合'];
wall.heatingHistoryMethod=history.method;
wall.coupling.externalStateSolves=stationSolves;
wall.coupling.boundaryEvaluations=boundaryEvaluations;
wall.coupling.configurationSource='calcFlightHeatingHistory 保存的气动热配置';
wall.flightHistory=struct('time',wall.time);
for key={'altitude','Mach','alpha','beta'}
    wall.flightHistory.(key{1})=interp1(flight.time,flight.(key{1}),wall.time,'linear');
end

%% 3. 驻点是独立材料点；仅在位置不随工况移动时提供单点壁温历程
meta=wall.coupling.boundaryMetadata;
point=struct('available',false,'status','驻点边界条件未全部可用');
pointValid=cellfun(@(m)~isempty(m) && m.stagnation.available,meta);
if all(pointValid)
    positions=cell2mat(cellfun(@(m)m.stagnation.position(:)',meta,'UniformOutput',false)');
    scale=max(1,max(max(mesh.V)-min(mesh.V)));
    fixedPosition=all(vecnorm(positions-positions(1,:),2,2)<=1e-8*scale);
    if fixedPosition
        point=stationaryPointHistory(meta,wall.time,initial,material,cfg,history.method);
    else
        point.status=['驻点位置随工况移动，不把不同材料位置拼成单点温度历程；' ...
            '请查看固定面元温度曲线'];
    end
end
wall.stagnationPoint=point;
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
mesh.wallThermal=wall;
% 初始壁温基准历程与实际热壁结果分别保存，避免把冷壁覆盖率当成热壁覆盖率。
mesh.flightHeatingHistory=history;

    function heat=evaluate(temperature,atTime)
        % 缓存仅保留最近一种飞行条件；恒定工况可跨步复用外流，
        % 不复用旧壁温热流。拒绝步回退到子步时会重建对应工况。
        condition=stationCfg.flight;
        for field={'altitude','Mach','alpha','beta'}
            condition.(field{1})=interp1(flight.time,flight.(field{1}),atTime,'linear');
        end
        values=[condition.altitude condition.Mach condition.alpha condition.beta];
        if isempty(snapshot) || ~isequal(values,cachedCondition)
            snapshot=calcFlightStationState(mesh,condition,stationCfg);
            flowSnapshot=captureAerothermalFlowSnapshot(snapshot);
            cachedCondition=values; currentCondition=condition;
            stationSolves=stationSolves+1;
        end
        localCfg=stationCfg; localCfg.flight=currentCondition;
        updated=recomputeSurfaceHeatingAtWallTemperature(snapshot, ...
            snapshot.externalFlow.freestream,localCfg,temperature);
        boundaryEvaluations=boundaryEvaluations+1;
        heat=updated.heating;
        heat.coverageReasonCode=updated.heatingCoverage.reasonCode;
        heat.boundaryMetadata=struct('time',atTime,'flight',condition, ...
            'flowSnapshot',flowSnapshot, ...
            'externalFlowMethod',snapshot.externalFlow.method, ...
            'axisymmetricMarchApplied',snapshot.axisymmetricMarch.applied, ...
            'axisymmetricMarchStatus',snapshot.axisymmetricMarch.status, ...
            'stagnation',snapshot.stagnationHeating, ...
            'thermalPropertyMethod',heat.thermalPropertyMethod, ...
            'computedAreaFraction',updated.heatingCoverage.summary.computedAreaFraction);
    end
end

function point=stationaryPointHistory(meta,time,initial,material,cfg,historyMethod)
% 单点仍用 Sutton-Graves 驻点边界及焓差反馈，不能冒充侧壁 Head 积分。
nt=numel(time);
h=cellfun(@(m)m.stagnation.heatTransferCoefficient,meta);
Taw=cellfun(@(m)m.stagnation.adiabaticWallTemperature,meta);
profile=struct('time',time,'initialWallTemperature',initial, ...
    'heatTransferCoefficient',h,'recoveryTemperature',Taw, ...
    'eligibleFace',true(1,nt));
modes=cellfun(@(m)m.thermalPropertyMethod,meta,'UniformOutput',false);
if all(strcmp(modes,'thermallyPerfect'))
    reference=Taw+0.5*(initial-Taw);
    mu=1.458e-6*reference.^(3/2)./(reference+110.4);
    profile.referenceLaw=struct('edgeTemperature',Taw,'edgePressure',ones(1,nt), ...
        'exponent',zeros(1,nt),'referenceViscosity0',mu.*287.05.*reference, ...
        'gasConstant',287.05,'thermalPropertyMode','thermallyPerfect', ...
        'referenceCp0',airEnthalpySecantCp(Taw,initial));
elseif any(strcmp(modes,'thermallyPerfect'))
    point=struct('available',false,'status','驻点热物性模型跨站改变，单点热历程未计算');
    return;
end
single=calcWallThermalHistory(struct('F',[1 2 3]),profile,material,cfg);
w=single.wallThermal;
point=struct('available',true,'status','固定位置驻点单点壁温；不是鼻部温度分布', ...
    'method',[historyMethod ' + ' w.method],'scope',w.scope, ...
    'time',time,'temperature',w.temperature(1,:), ...
    'convectiveHeatFlux',w.convectiveHeatFlux(1,:), ...
    'radiativeHeatFlux',w.radiativeHeatFlux(1,:), ...
    'netHeatFlux',w.netHeatFlux(1,:),'oneDimensional',logical(w.oneDimensionalFace), ...
    'backTemperature',w.backTemperature(1,:),'backHeatFlux',w.backHeatFlux(1,:), ...
    'backConvectiveHeatFlux',w.backConvectiveHeatFlux(1,:), ...
    'backRadiativeHeatFlux',w.backRadiativeHeatFlux(1,:));
end
