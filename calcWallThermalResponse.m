function mesh = calcWallThermalResponse(mesh,timeSeconds,material,cfg)
%CALCWALLTHERMALRESPONSE 固定工况的壁温-边界层耦合及低阶对照模型。
% 有完整外流和准流线时默认每个隐式步重算沿程边界层，所有有效
% 面元求一维厚度导热；cfg.heating.wallCoupling='local' 可用旧对照。
% 下述 Biot 分流和幂指数反馈仅属于 local 分支或给定 h/单点模型。
% 能量平衡（单位面积）：rho_m*c_m*L*dTw/dt =
%   h*(Taw-Tw) - emissivity*sigma*(Tw^4-Tr^4)。
% 每一步用隐式 Euler + Newton 迭代，处理 T^4 辐射非线性。
% Bi=h*(L/2)/k_m <=0.1 使用集中热容；更大的 Bi 自动转入
% calcWallConduction1D 的一维厚度方向有限体积模型。
% 有内置外流与热流分区时，按当前壁温更新 Eckert 参考性质及换热系数；
% 温度相关空气模式还更新焓差等效比热及 Newton 导数。
% 外流、恢复温度、转捩类别和压力梯度仍固定在给定飞行工况。
% 材料参数由调用方提供，不从 STL 自动猜测。

%% 1. 检查时间、材料和已计算的换热边界条件
if nargin<4, cfg=struct; end
% 单工况剖面的换热割线不能直接作为旧热历史反馈规律使用。
if isfield(mesh,'heating') && isfield(mesh.heating,'historyCompatible') && ...
        ~mesh.heating.historyCompatible
    error('RocketAeroHeat:WallThermal:SingleConditionModel', ...
        '当前球鼻局部相似模型仅验证单工况热流；壁体热历史耦合尚未验证。');
end
backBoundary=wallBackBoundary(material);
if ~isfield(mesh,'heating') || ~isfield(mesh.heating,'heatTransferCoefficient') || ...
        ~isfield(mesh.heating,'recoveryTemperature')
    error('RocketAeroHeat:WallThermal:Heating','请先计算包含换热系数的面元热流。');
end
validateattributes(timeSeconds,{'numeric'}, ...
    {'real','finite','vector','nonempty','nonnegative'},mfilename,'timeSeconds');
t=double(timeSeconds(:)');
if t(1)~=0 || any(diff(t)<=0)
    error('RocketAeroHeat:WallThermal:Time','时间须从 0 开始并严格递增。');
end
layeredMaterial=isfield(material,'layers') && ~isempty(material.layers);
if layeredMaterial
    if ~isfield(cfg,'flight') || ~isfield(cfg.flight,'Twall')
        error('RocketAeroHeat:LayeredWall:Initial','多层材料需要 cfg.flight.Twall 初温。');
    end
    material=normalizeLayeredWallMaterial(material,cfg.flight.Twall);
end
for name={'density','specificHeat','thickness','thermalConductivity', ...
        'emissivity','radiationTemperature'}
    key=name{1};
    if ~isfield(material,key)
        error('RocketAeroHeat:WallThermal:Material','缺少 material.%s。',key);
    end
    validateattributes(material.(key),{'numeric'}, ...
        {'real','finite','scalar'},mfilename,['material.' key]);
end
for name={'density','specificHeat','thickness','thermalConductivity', ...
        'radiationTemperature'}
    if material.(name{1})<=0
        error('RocketAeroHeat:WallThermal:Material','%s 必须为正。',name{1});
    end
end
if material.emissivity<0 || material.emissivity>1
    error('RocketAeroHeat:WallThermal:Material','发射率应在 0~1。');
end
if ~isfield(cfg,'flight') || ~isfield(cfg.flight,'Twall')
    error('RocketAeroHeat:WallThermal:Initial','需要 cfg.flight.Twall 初始壁温。');
end
initial=cfg.flight.Twall;
validateattributes(initial,{'numeric'},{'real','finite','scalar','positive'});
materialProperties=wallMaterialProperties(material,initial);
material.specificHeat=materialProperties.specificHeat;
material.thermalConductivity=materialProperties.thermalConductivity;
variableMaterial=isfield(material,'propertyTable') && ~isempty(material.propertyTable);
if isfield(mesh.heating,'wallTemperatureByFace') && ...
        any(~isfinite(mesh.heating.wallTemperatureByFace) | ...
        abs(mesh.heating.wallTemperatureByFace-initial)>1e-9*max(1,initial))
    error('RocketAeroHeat:WallThermal:StaleHeating', ...
        '当前热流使用非均匀或未知壁温，不能当作统一初始壁温的热边界。');
end
if isfield(mesh.heating,'wallTemperature') && ...
        abs(mesh.heating.wallTemperature-initial)>1e-9*max(1,initial)
    error('RocketAeroHeat:WallThermal:StaleHeating', ...
        '热流边界条件使用的初始壁温与当前 cfg.flight.Twall 不一致。');
end
h=mesh.heating.heatTransferCoefficient(:);
Taw=mesh.heating.recoveryTemperature(:);
nf=size(mesh.F,1); nt=numel(t);
if numel(h)~=nf || numel(Taw)~=nf
    error('RocketAeroHeat:WallThermal:Size','换热边界条件与面元数不匹配。');
end
if nf*nt>1e7
    error('RocketAeroHeat:WallThermal:Memory','输出矩阵过大，请减少时间采样点。');
end

%% 2. 从当前外流建立壁温相关的局部参考性质反馈
% 完整内置外流路径默认进行沿程耦合；简化给定 h/单点模型保留局部反馈。
couplingMode='streamline';
if isfield(cfg,'heating') && isfield(cfg.heating,'wallCoupling')
    couplingMode=char(cfg.heating.wallCoupling);
end
if ~ismember(couplingMode,{'streamline','local'})
    error('RocketAeroHeat:WallThermal:Coupling','wallCoupling 须为 streamline 或 local。');
end
useCoupled=strcmp(couplingMode,'streamline') && ...
    isfield(mesh,'externalFlow') && isfield(mesh,'surfaceStreamlines');
if useCoupled
    evaluator=@(temperature)evaluateCoupledBoundary(mesh,cfg,temperature);
    initialBoundary=evaluator(initial*ones(nf,1));
    result=solveCoupledWallConduction(evaluator,initialBoundary,t,material,initial);
else
law=[];
if isfield(mesh,'boundaryLayerEdge') && ...
        isfield(mesh.heating,'laminarFace') && ...
        isfield(mesh.heating,'turbulentFace') && ...
        isfield(mesh.heating,'Prandtl')
    e=mesh.boundaryLayerEdge;
    if isfield(e,'temperature') && isfield(e,'pressure') && ...
            isfield(e,'freestream') && isfield(e.freestream,'R')
        exponent=nan(nf,1);
        exponent(mesh.heating.laminarFace)=0.5;
        exponent(mesh.heating.turbulentFace)=0.2;
        Te=e.temperature(:); p=e.pressure(:);
        % 鼻端延拓沿用驻点固定 h；没有逐面元外缘状态，指数取 0。
        noseFace=false(nf,1);
        if isfield(mesh.heating,'noseFace')
            noseFace=mesh.heating.noseFace(:);
            exponent(noseFace)=0;
            Te(noseFace)=Taw(noseFace);
            p(noseFace)=e.freestream.pressure;
        end
        reference=Te+0.5*(initial-Te)+0.22*(Taw-Te);
        mu=1.458e-6*reference.^(3/2)./(reference+110.4);
        nu0=mu.*e.freestream.R.*reference./p;
        eligible=mesh.heating.eligibleFace(:);
        if any(~isfinite(nu0(eligible)) | nu0(eligible)<=0 | ...
                ~isfinite(exponent(eligible)))
            error('RocketAeroHeat:WallThermal:Reference', ...
                '外流参考性质或转捩类别缺失，无法更新换热系数。');
        end
        law=struct('edgeTemperature',Te,'edgePressure',p, ...
            'exponent',exponent,'referenceViscosity0',nu0, ...
            'gasConstant',e.freestream.R);
        if isfield(mesh.heating,'thermalPropertyMethod') && ...
                strcmp(mesh.heating.thermalPropertyMethod,'thermallyPerfect')
            law.thermalPropertyMode='thermallyPerfect';
            law.referenceCp0=airEnthalpySecantCp(Taw,initial);
            if any(~isfinite(law.referenceCp0(eligible)) | ...
                    law.referenceCp0(eligible)<=0)
                error('RocketAeroHeat:WallThermal:AirRange', ...
                    '初始壁温或恢复温度超出干空气焓模型适用域。');
            end
        end
    end
end

%% 3. 检查集中热容 Biot 数；按壁温范围估计最大的 Biot 数
Bi=h*(material.thickness/2)/material.thermalConductivity;
hUpper=h;
if ~isempty(law)
    validEdge=mesh.heating.eligibleFace(:) & isfinite(h) & isfinite(Taw);
    upperTemperature=max(initial,Taw(validEdge));
    hUpper(validEdge)=evaluateWallConvection(h(validEdge),Taw(validEdge), ...
        upperTemperature,subsetLaw(law,validEdge));
end
BiUpper=max(Bi,hUpper*(material.thickness/2)/material.thermalConductivity);
valid=mesh.heating.eligibleFace(:) & isfinite(h) & h>=0 & ...
    isfinite(Taw) & Taw>0 & BiUpper<=0.1;
conduction=mesh.heating.eligibleFace(:) & isfinite(h) & h>=0 & ...
    isfinite(Taw) & Taw>0 & BiUpper>0.1;
if variableMaterial
    % Biot 数不再决定变物性材料是否简化为常比热集中热容。
    conduction=valid|conduction; valid(:)=false;
    BiUpper=hUpper*(material.thickness/2)/min(material.propertyTable.thermalConductivity);
end
if layeredMaterial || backBoundary.active
    conduction=mesh.heating.eligibleFace(:) & isfinite(h) & h>=0 & isfinite(Taw) & Taw>0;
    valid(:)=false; Bi(:)=NaN; BiUpper(:)=NaN;
end
C=material.density*material.specificHeat*material.thickness;
sigma=5.670374419e-8;
epsWall=material.emissivity;
Tr=material.radiationTemperature;
temperature=nan(nf,nt);
qConv=nan(nf,nt); qRad=qConv; qNet=qConv; hHistory=qConv;
temperature(valid,1)=initial;
hHistory(valid,1)=h(valid);
qConv(valid,1)=h(valid).*(Taw(valid)-initial);
qRad(valid,1)=epsWall*sigma*(initial^4-Tr^4);
qNet(valid,1)=qConv(valid,1)-qRad(valid,1);

%% 4. 隐式时间步进；每次 Newton 迭代都用当前壁温更新 h
ids=find(valid);
lumpedLaw=subsetLaw(law,valid);
for j=2:nt
    dt=t(j)-t(j-1);
    old=temperature(ids,j-1);
    value=old;
    for iteration=1:30
        [currentH,derivativeH]=evaluateWallConvection( ...
            h(ids),Taw(ids),value,lumpedLaw);
        residual=C/dt*(value-old)-currentH.*(Taw(ids)-value)+ ...
            epsWall*sigma*(value.^4-Tr^4);
        derivative=C/dt+currentH-derivativeH.*(Taw(ids)-value)+ ...
            4*epsWall*sigma*value.^3;
        correction=residual./derivative;
        value=max(1,value-correction);
        if all(abs(correction)<=1e-9*max(1,value)), break; end
    end
    if any(~isfinite(value)) || any(abs(correction)>1e-6*max(1,value))
        error('RocketAeroHeat:WallThermal:Convergence', ...
            '隐式壁温迭代未收敛，请检查时间步长和材料输入。');
    end
    temperature(ids,j)=value;
    currentH=evaluateWallConvection( ...
        h(ids),Taw(ids),value,lumpedLaw);
    hHistory(ids,j)=currentH;
    qConv(ids,j)=currentH.*(Taw(ids)-value);
    qRad(ids,j)=epsWall*sigma*(value.^4-Tr^4);
    qNet(ids,j)=qConv(ids,j)-qRad(ids,j);
end

%% 5. 对集中热容不适用的面元求一维厚度导热
oneDimensional=struct;
if any(conduction)
    conductionLaw=subsetLaw(law,conduction);
    oneDimensional=calcWallConduction1D( ...
        h(conduction),Taw(conduction),t,material,initial, ...
        conductionLaw);
    temperature(conduction,:)=oneDimensional.surfaceTemperature;
    for j=1:nt
        currentH=evaluateWallConvection(h(conduction),Taw(conduction), ...
            temperature(conduction,j),conductionLaw);
        hHistory(conduction,j)=currentH;
        qConv(conduction,j)=currentH.* ...
            (Taw(conduction)-temperature(conduction,j));
        qRad(conduction,j)=epsWall*sigma* ...
            (temperature(conduction,j).^4-Tr^4);
        qNet(conduction,j)=qConv(conduction,j)-qRad(conduction,j);
    end
end

%% 6. 记录算法、材料、适用性和完整时间历程
result=struct;
if isempty(law)
    result.method=['Bi<=0.1：集中热容；Bi>0.1：一维有限体积导热；' ...
        '隐式 Euler 与表面辐射迭代，给定常数 h'];
else
    result.method=['Eckert 参考性质壁温反馈 + Biot 分流 + 隐式 Euler；' ...
        '厚壁采用一维有限体积导热'];
end
result.scope=['均匀常物性材料、绝热背面；外流、转捩和恢复温度固定，' ...
    '侧壁局部参考黏度随壁温更新，层流/湍流指数 1/2 与 1/5 为局部近似，' ...
    '不随壁温重积分 Head 状态或重判分离；鼻端延拓保持固定 h；不含烧蚀'];
result.time=t; result.material=material;
result.initialTemperature=initial;
result.temperature=temperature;
result.convectiveHeatFlux=qConv;
result.radiativeHeatFlux=qRad;
result.netHeatFlux=qNet;
result.heatTransferCoefficientHistory=hHistory;
result.Biot=Bi; result.maximumEstimatedBiot=BiUpper;
result.wallTemperatureFeedback=~isempty(law);
result.validFace=valid|conduction;
result.lumpedFace=valid; result.oneDimensionalFace=conduction;
result.oneDimensional=oneDimensional;
result.materialMethod=materialProperties.method;
if variableMaterial
    result.method=[result.method '；变物性材料全部采用 cp(T) 积分储热/Kirchhoff 一维导热'];
    result.scope=strrep(result.scope,'均匀常物性材料','单层常密度、用户温度相关 cp/k 材料');
end
if layeredMaterial
    result.method='逐层有限体积 + cp(T) 积分储热/Kirchhoff 通量 + 接触热阻 + 隐式 Euler/辐射';
    result.materialMethod='逐层独立材料物性与界面热阻；不使用平均物性求解';
    result.scope=strrep(result.scope,'均匀常物性材料','多层常密度材料、各层可给 cp(T)/k(T)');
end
result.summary=struct('calculatedFaces',nnz(valid|conduction), ...
    'lumpedFaces',nnz(valid),'oneDimensionalFaces',nnz(conduction));
if any(valid|conduction)
    result.summary.maximumFinalTemperature=max(temperature(valid|conduction,end));
else
    result.summary.maximumFinalTemperature=NaN;
end

%% 7. 球形钝鼻驻点单独求一个材料热响应；不向邻近面元外推
end % 局部反馈分支；以下驻点单点记录由两种模型共用。
% 耦合求解可能自动插入子步；驻点和面元历史必须共享实际接受的时间轴。
t=result.time;
result=completeWallEnergyBalance(result);
point=struct('available',false,'status','当前工况没有可用的驻点单点热流');
if isfield(mesh,'stagnationHeating') && mesh.stagnationHeating.available
    stag=mesh.stagnationHeating;
    if ~isfield(stag,'wallTemperature') || ...
            ~isfinite(stag.wallTemperature) || ...
            abs(stag.wallTemperature-initial)>1e-9*max(1,initial)
        error('RocketAeroHeat:WallThermal:StaleStagnation', ...
            '驻点热流使用的初始壁温与当前工况不一致。');
    end
    if isfinite(stag.heatTransferCoefficient) && ...
            stag.heatTransferCoefficient>0 && ...
            isfinite(stag.adiabaticWallTemperature) && ...
            stag.adiabaticWallTemperature>0
        single=struct('F',[1 2 3]);
        single.heating=struct('heatTransferCoefficient', ...
            stag.heatTransferCoefficient, ...
            'recoveryTemperature',stag.adiabaticWallTemperature, ...
            'eligibleFace',true,'wallTemperature',initial);
        if isfield(stag,'thermalPropertyMode') && ...
                strcmp(stag.thermalPropertyMode,'thermallyPerfect')
            % 驻点单点的空气焓差随壁温更新；黏度指数置零。
            single.heating.thermalPropertyMethod='thermallyPerfect';
            single.heating.laminarFace=false;
            single.heating.turbulentFace=false;
            single.heating.noseFace=true;
            single.heating.Prandtl=0.72;
            % 指数为零时参考压力和气体常数在黏度比中严格相消。
            single.boundaryLayerEdge=struct('temperature', ...
                stag.adiabaticWallTemperature,'pressure',1, ...
                'freestream',struct('R',287.05,'pressure',1));
        end
        single=calcWallThermalResponse(single,t,material,cfg);
        local=single.wallThermal;
        point=struct('available',true, ...
            'status','单个驻点的材料热响应；不是鼻部温度分布', ...
            'method',[stag.method ' + ' local.method], ...
            'scope',[stag.scope '；' local.scope], ...
            'position',stag.position,'radius',stag.radius, ...
            'time',t,'temperature',local.temperature(1,:), ...
            'convectiveHeatFlux',local.convectiveHeatFlux(1,:), ...
            'radiativeHeatFlux',local.radiativeHeatFlux(1,:), ...
            'netHeatFlux',local.netHeatFlux(1,:), ...
            'heatTransferCoefficient',stag.heatTransferCoefficient, ...
            'oneDimensional',logical(local.oneDimensionalFace), ...
            'backTemperature',local.backTemperature(1,:), ...
            'backHeatFlux',local.backHeatFlux(1,:), ...
            'backConvectiveHeatFlux',local.backConvectiveHeatFlux(1,:), ...
            'backRadiativeHeatFlux',local.backRadiativeHeatFlux(1,:));
    end
end

result.stagnationPoint=point;
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
mesh.wallThermal=result;
end

function part=subsetLaw(law,mask)
% 仅截取当前求解面元；空 law 表示保持用户给定的固定 h。
if isempty(law), part=[]; return; end
part=law;
for name={'edgeTemperature','edgePressure','exponent','referenceViscosity0'}
    values=law.(name{1});
    part.(name{1})=values(mask);
end
if isfield(law,'referenceCp0'), part.referenceCp0=law.referenceCp0(mask); end
end

function heat=evaluateCoupledBoundary(mesh,cfg,temperature)
% 用当前面元壁温重算整条热边界，附带覆盖原因以审计途中失效。
updated=recomputeSurfaceHeatingAtWallTemperature(mesh, ...
    mesh.externalFlow.freestream,cfg,temperature);
heat=updated.heating;
heat.coverageReasonCode=updated.heatingCoverage.reasonCode;
end
