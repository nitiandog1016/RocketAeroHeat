function result=solveSphericalNoseLaminarHeating(mesh,freestream,wallTemperature,cfg)
%SOLVESPHERICALNOSELAMINARHEATING 已核验球鼻的连续流层流区域热边界。
% 模块作用：从同一解析球面弧线上的压力/驻点熵外缘积分xi，逐站解
% Levy-Lees局部相似速度/总焓剖面，直接求热流；不采用驻点值*cos(theta)。
% 恢复温度由同一剖面的q(Twall)=0求根；换热系数为q/(Taw-Twall)割线值。
% 来源：NASA SP-8062 §2.3.2、4.2.1；NASA TP-3271式6、18–25、56。
% 本版仍忽略沿xi的非相似项，不解熵层吞入/横向涡量、转捩或区域接口。
% 只在多尺度球面核验的完整光滑前缀上计算；后体不由该解外推补齐。
% 仅接受实际来流沿-X的零合成攻角；外缘模块的斜来流几何范围并不
% 构成本轴对称热边界的适用依据，不以轴对称剖面代替斜流横向解。
% maxPolarAngle由现有球鼻外缘配置筛选，Re与delta/R是显式项目筛选，
% 不是来源所给的普适准确度界限。NASA7越域不回退定比热。
if nargin<4, cfg=struct; end
nf=size(mesh.F,1);
validateattributes(wallTemperature,{'numeric'},{'real','vector'});
Tw=double(wallTemperature(:));
if isscalar(Tw), Tw=repmat(Tw,nf,1); end
if numel(Tw)~=nf || any(isinf(Tw))
    error('RocketAeroHeat:SphereHeating:Wall','壁温须为标量或逐面元向量；未知边界可保留NaN。');
end
edge=solveSphericalNoseEdge(mesh,freestream,cfg);
nanField=nan(nf,1);
result=struct('available',false,'faceMask',false(nf,1), ...
    'heatFlux',nanField,'heatTransferCoefficient',nanField, ...
    'recoveryTemperature',nanField,'wallTemperatureInput',Tw, ...
    'wallShearStress',nanField,'massThickness',nanField,'delta',nanField, ...
    'xi',nanField,'beta',nanField,'pressureGradient',nanField, ...
    'velocityGradient',nanField,'inputMomentumRelativeResidual',nanField, ...
    'adiabaticHeatResidual',nanField,'polarAngleDeg',edge.polarAngleDeg, ...
    'statusCode',ones(nf,1,'uint8'),'edgeState',edge, ...
    'statusNames',{{'已解球鼻局部相似工程热边界','几何/外缘不适用', ...
        '连续介质筛选未通过','壁温未知或物性域外','层流Re筛选未通过', ...
        '剖面或绝热热边界未取得有效解','delta/R薄层筛选未通过', ...
        '非零合成攻角不在本轴对称热边界范围'}}, ...
    'method','球鼻驻点熵外缘 + Levy-Lees局部相似BVP壁面导数', ...
    'scope',['零合成攻角、已核验球面前缀的连续流附着层流工程近似；' ...
        '不解沿程非相似、熵层吞入、横向涡量、交界波系或实验精度'], ...
    'heatFluxMethod','NASA TP-3271式56，同一BVP总焓壁面导数', ...
    'recoveryTemperatureMethod','同一BVP的q(Twall)=0绝热热边界求根，非给定总温', ...
    'heatTransferCoefficientMethod','同一求解热流与绝热壁温的割线q/(Taw-Twall)', ...
    'sourceURLs',{{'https://ntrs.nasa.gov/api/citations/19710021703/downloads/19710021703.pdf', ...
        'https://ntrs.nasa.gov/api/citations/19930006245/downloads/19930006245.pdf'}}, ...
    'nonSimilarTermsSolved',false,'entropySwallowingSolved',false, ...
    'interfaceSolved',false,'experimentalValidation',false,'historyCompatible',false, ...
    'noseRadius',edge.noseRadius,'angularStepDeg',2,'maximumDistinctWallTemperatures',16, ...
    'maximumDeltaRadiusFraction',.05,'transitionRe',5e5,'maximumContinuumKnudsen',.01, ...
    'freestreamKnudsen',NaN,'profileCalls',0,'uniqueStations',0,'stationResults',{{}}, ...
    'incidenceDeg',atan2d(norm(freestream.unitFlowDirection(2:3)), ...
        -freestream.unitFlowDirection(1)), 'maximumIncidenceDeg',1e-8, ...
    'status','尚无有效球鼻区域热边界');
if ~edge.available, result.status=edge.status; return; end
if result.incidenceDeg>result.maximumIncidenceDeg
    result.statusCode(edge.geometryFace)=7;
    result.status='非零合成攻角：本版轴对称局部相似热边界不适用；未解横流或区域接口'; return;
end
R=edge.noseRadius;
if isfield(cfg,'heating') && isfield(cfg.heating,'transitionRe')
    result.transitionRe=cfg.heating.transitionRe;
end
if isfield(cfg,'heating') && isfield(cfg.heating,'maxContinuumKnudsen')
    result.maximumContinuumKnudsen=cfg.heating.maxContinuumKnudsen;
end
if isfield(cfg,'heating') && isfield(cfg.heating,'sphericalNoseAngleStepDeg')
    result.angularStepDeg=cfg.heating.sphericalNoseAngleStepDeg;
end
validateattributes(result.angularStepDeg,{'numeric'},{'scalar','real','finite','>=',.25,'<=',5});
validateattributes(result.transitionRe,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(result.maximumContinuumKnudsen,{'numeric'}, ...
    {'scalar','real','finite','>',0,'<=',.01});
muInf=mu(freestream.temperature);
result.freestreamKnudsen=airMeanFreePath(freestream.pressure, ...
    freestream.temperature,muInf,freestream.R)/R;
if ~isfinite(result.freestreamKnudsen) || result.freestreamKnudsen>=result.maximumContinuumKnudsen
    result.statusCode(edge.algebraicallyValid)=2;
    result.status='连续介质筛选未通过；未计算球鼻无滑移层流热边界'; return;
end

%% 1. 同一球面压力及正激波熵状态；导数和驻点二阶曲率互相兼容
normal=edge.normalShock; p0=normal.stagnationPressure; T0=normal.stagnationTemperature;
rho0=p0/(freestream.R*T0); pressureDrop=p0-freestream.pressure;
a=sqrt(2*pressureDrop/(rho0*R^2));
result.stagnationVelocityGradient=a;
result.stagnationPressureCurvature=-2*pressureDrop/R^2;
result.stagnationMomentumCurvatureResidual= ...
    abs(result.stagnationPressureCurvature+rho0*a^2)/(rho0*a^2);
result.stagnationPressureCurvatureVerified=result.stagnationMomentumCurvatureResidual<1e-12;
mode=edge.thermalPropertyMode;
if strcmp(mode,'thermallyPerfect'), [~,H0]=airThermalProperties(T0);
else, cp=freestream.gamma*freestream.R/(freestream.gamma-1); H0=cp*T0; end
in=edge.algebraicallyValid & isfinite(Tw) & Tw>=80 & Tw<=2500;
result.statusCode(edge.algebraicallyValid & ~in)=3;
ids=find(in);
if isempty(ids), result.status='没有有效壁温与球鼻外缘相匹配的面元'; return; end
walls=unique(Tw(ids));
if numel(walls)>result.maximumDistinctWallTemperatures
    result.statusCode(ids)=5;
    result.status='壁温分组超过16；本版只支持少量给定壁温的单时刻角度表，不用粗壁温插值补齐'; return;
end
% 每个给定壁温只建一张角度表；数万面元不会变成数万个BVP。
% 只在两个有效相邻节点之间线性插值，不跨过未解节点或几何区域。
maxAngle=max(edge.polarAngleDeg(ids));
angles=unique([0:result.angularStepDeg:maxAngle,maxAngle])';
repeatedWalls=repelem(walls(:),numel(angles));
uniqueKeys=[repmat(angles,numel(walls),1),repeatedWalls(:)];
result.uniqueStations=size(uniqueKeys,1);
result.stationResults=cell(result.uniqueStations,1);
tableFields={'heatFlux','heatTransferCoefficient','recoveryTemperature','wallShearStress', ...
    'massThickness','delta','xi','beta','pressureGradient','velocityGradient', ...
    'inputMomentumRelativeResidual','adiabaticHeatResidual'};
nodeValues=nan(result.uniqueStations,numel(tableFields));
nodeCode=5*ones(result.uniqueStations,1,'uint8');
result.profileNodeIDs=nan(nf,2); result.profileNodeWeight=nanField;
result.angularInterpolation='仅在同壁温的两个有效相邻BVP角度节点间线性插值；不跨未解节点';
result.heatFluxMethod=[result.heatFluxMethod '；' result.angularInterpolation];
result.recoveryTemperatureMethod=[result.recoveryTemperatureMethod '；' result.angularInterpolation];

%% 2. 从驻点沿真实球面弧线积分xi，逐站求壁面热流和绝热边界
for j=1:result.uniqueStations
    phi=deg2rad(uniqueKeys(j,1)); wall=uniqueKeys(j,2);
    state=stateAt(phi); s=R*phi; h=R*sin(phi);
    if ~state.valid
        nodeCode(j)=5; continue;
    end
    Re=state.density*state.speed*s/mu(state.temperature);
    if Re>=result.transitionRe
        nodeCode(j)=4; continue;
    end
    if phi==0
        xi=0; beta=.5; dP=0; dU=a;
    else
        % 变量替换t=phi*z，避免很小鼻端站的xi~phi^4被绝对积分容差淹没。
        xi=R*phi*integral(@(z)xiIntegrand(phi*z),0,1, ...
            'RelTol',1e-8,'AbsTol',0);
        dP=-pressureDrop*sin(2*phi)/R;
        dU=-dP/(state.density*state.speed);
        beta=2*xi*dU/(state.density*mu(state.temperature)*state.speed^2*h^2);
    end
    supplied=struct('temperature',state.temperature,'pressure',state.pressure, ...
        'density',state.density,'speed',state.speed,'totalEnthalpy',H0,'sources',struct);
    for key={'temperature','density','speed','totalEnthalpy'}
        supplied.sources.(key{1})='同一正激波驻点熵和总焓，在球面修正牛顿静压等熵恢复';
    end
    supplied.sources.pressure='解析球面p=pinf+(p02-pinf)cos(phi)^2修正牛顿工程压力';
    station=struct('streamlineMetric',h,'xi',xi,'beta',beta, ...
        'pressureGradient',dP,'velocityGradient',dU,'sources',struct);
    station.sources.streamlineMetric='解析球面流管度量h=Rn*sin(phi)';
    station.sources.xi='NASA TP-3271式6，驻点起算积分rho_e*mu_e*U_e*h^2 ds';
    station.sources.beta='NASA TP-3271式23，由同一xi及dUe/ds计算';
    station.sources.pressureGradient='同一解析球面修正牛顿p对s=Rn*phi求导';
    station.sources.velocityGradient='固定驻点熵外缘：dUe/ds=-pe_s/(rhoe*Ue)；驻点取解析曲率极限';
    station.sources.wallTemperature='给定逐面元壁温边界，不是本模块求解壁温';
    if phi==0
        station.stagnationVelocityGradient=a;
        station.sources.stagnationVelocityGradient='由同一球面压力二阶曲率得a=sqrt(2*(p02-pinf)/(rhos*Rn^2))';
    end
    current=solveAtWall(wall);
    result.stationResults{j}=struct('polarAngleDeg',rad2deg(phi),'wallTemperature',wall, ...
        'wallProfile',current,'edge',supplied,'station',station);
    if ~current.available
        nodeCode(j)=5; continue;
    end
    if current.delta/R>result.maximumDeltaRadiusFraction
        nodeCode(j)=6; continue;
    end
    try
        if phi==0
            taw=T0; adiabatic=solveAtWall(taw);
            bracket=[T0 T0]; trials=[];
        else
            [bracket,trials]=adiabaticBracket();
            taw=fzero(@heatAtWall,bracket,optimset('TolX',1e-7,'Display','off'));
            adiabatic=solveAtWall(taw);
        end
        if ~adiabatic.available || abs(adiabatic.heatFlux)>1e-5*max(abs(current.heatFlux),1)
            error('RocketAeroHeat:SphereHeating:Adiabatic','绝热壁温热流残差检查未通过。');
        end
        if abs(taw-wall)>1e-5
            coefficient=current.heatFlux/(taw-wall);
        else
            % 已在绝热壁温附近；用可审查的1 K热边界导数极限求h，
            % 不由0/0补零，不裁剪真实q。1 K为显式数值差分尺度。
            lower=max(80,taw-1); upper=min(2500,taw+1);
            low=solveAtWall(lower); high=solveAtWall(upper);
            if ~low.available || ~high.available, error('RocketAeroHeat:SphereHeating:Adiabatic','热导数未解。'); end
            coefficient=-(high.heatFlux-low.heatFlux)/(upper-lower);
        end
        if ~isfinite(coefficient) || coefficient<=0
            error('RocketAeroHeat:SphereHeating:Adiabatic','局部热边界系数不是有限正值。');
        end
    catch exception
        result.stationResults{j}.thermalMessage=exception.message;
        nodeCode(j)=5; continue;
    end
    result.stationResults{j}.adiabaticProfile=adiabatic;
    result.stationResults{j}.adiabaticBracket=bracket;
    result.stationResults{j}.adiabaticBracketTrials=trials;
    result.stationResults{j}.recoveryTemperature=taw;
    nodeValues(j,:)=[current.heatFlux,coefficient,taw,current.wallShearStress, ...
        current.massThickness,current.delta,xi,beta,dP,dU, ...
        current.inputMomentumRelativeResidual,adiabatic.heatFlux];
    nodeCode(j)=0;
end
%% 3. 将区域热边界映射到真实面元；每一插值仍可追溯两个BVP节点
for k=1:numel(ids)
    id=ids(k); candidates=find(uniqueKeys(:,2)==Tw(id));
    lo=find(uniqueKeys(candidates,1)<=edge.polarAngleDeg(id),1,'last');
    hi=find(uniqueKeys(candidates,1)>=edge.polarAngleDeg(id),1,'first');
    if isempty(lo) || isempty(hi), result.statusCode(id)=5; continue; end
    low=candidates(lo); high=candidates(hi);
    if nodeCode(low)~=0 || nodeCode(high)~=0
        result.statusCode(id)=max(nodeCode(low),nodeCode(high)); continue;
    end
    if low==high, weight=0;
    else, weight=(edge.polarAngleDeg(id)-uniqueKeys(low,1))/(uniqueKeys(high,1)-uniqueKeys(low,1)); end
    values=(1-weight)*nodeValues(low,:)+weight*nodeValues(high,:);
    for f=1:numel(tableFields), result.(tableFields{f})(id)=values(f); end
    if abs(result.recoveryTemperature(id)-Tw(id))>1e-5
        result.heatTransferCoefficient(id)=result.heatFlux(id)/(result.recoveryTemperature(id)-Tw(id));
    end
    result.profileNodeIDs(id,:)=[low high]; result.profileNodeWeight(id)=weight;
    result.faceMask(id)=true; result.statusCode(id)=0;
end
result.angularProfileTable=struct('polarAngleDeg',uniqueKeys(:,1), ...
    'wallTemperature',uniqueKeys(:,2),'statusCode',nodeCode,'fieldNames',{tableFields},'values',nodeValues);
result.available=any(result.faceMask);
result.status=sprintf('同球面弧线/BVP已计算%d面元、%d个独立局部站；区域接口仍未求解', ...
    nnz(result.faceMask),result.uniqueStations);
result.fieldMetadata=struct( ...
    'heatFlux',struct('valueKind','computed','method',result.heatFluxMethod,'scope',result.scope), ...
    'recoveryTemperature',struct('valueKind','computed','method',result.recoveryTemperatureMethod,'scope',result.scope), ...
    'wallTemperatureInput',struct('valueKind','prescribed','method','给定壁温热边界','scope','不作为求解壁温'), ...
    'heatTransferCoefficient',struct('valueKind','computed','method',result.heatTransferCoefficientMethod,'scope',result.scope));

    function r=solveAtWall(t)
        result.profileCalls=result.profileCalls+1;
        r=solveAxisymmetricLaminarMassProfile(freestream,supplied,t,station,cfg);
    end
    function q=heatAtWall(t)
        r=solveAtWall(t);
        if ~r.available, error('RocketAeroHeat:SphereHeating:Adiabatic','绝热壁温试探剖面无解。'); end
        q=r.heatFlux;
    end
    function [bracket,trials]=adiabaticBracket()
        % 模块作用：寻找同一有效附着剖面的绝热根括区。
        % Tw=T0可能产生非单调速度而被剖面内核拒绝；这不等于较低
        % 绝热壁温无解。仅缩小数值试探区间，寻找有效的异号端点，
        % 不放宽剖面门限，不使用无解热流，也不跨过拒绝节点插值。
        cold=solveAtWall(state.temperature);
        trials=[state.temperature,double(cold.available),double(cold.statusCode),cold.heatFlux];
        result.stationResults{j}.adiabaticBracketTrials=trials;
        if ~cold.available || cold.heatFlux<0
            error('RocketAeroHeat:SphereHeating:Adiabatic','无法建立有效的绝热根冷侧附着剖面。');
        end
        lower=state.temperature; upper=T0;
        probe=upper;
        for trial=1:24
            hot=solveAtWall(probe);
            trials(end+1,:)=[probe,double(hot.available),double(hot.statusCode),hot.heatFlux]; %#ok<AGROW>
            result.stationResults{j}.adiabaticBracketTrials=trials;
            if hot.available
                if hot.heatFlux<=0
                    bracket=[lower probe]; return;
                end
                if probe==T0
                    error('RocketAeroHeat:SphereHeating:Adiabatic','总温端点有效但热流仍为正，物理温度区间内未取得绝热根。');
                end
                lower=probe;
            else
                upper=probe;
            end
            probe=(lower+upper)/2;
        end
        error('RocketAeroHeat:SphereHeating:Adiabatic','未找到两端均有效且异号的绝热根区间；保留未解。');
    end
    function z=xiIntegrand(ph)
        e=stateAt(ph);
        if any(~e.valid), error('RocketAeroHeat:SphereHeating:Edge','xi积分中的外缘状态未解。'); end
        z=e.density.*mu(e.temperature).*e.speed.*(R*sin(ph)).^2;
    end
    function st=stateAt(ph)
        p=freestream.pressure+pressureDrop*cos(ph).^2;
        if strcmp(mode,'thermallyPerfect')
            % 既有恢复内核以列向量作热力反解；积分器传入行节点时，
            % 先以列调用再恢复节点形状，避免行/列隐式扩展成N×N。
            st=recoverThermallyPerfectEdgeState(p(:),repmat(p0,numel(p),1),T0,freestream.R);
            fields=fieldnames(st);
            for ii=1:numel(fields)
                key=fields{ii}; value=st.(key);
                if (isnumeric(value)||islogical(value)) && numel(value)==numel(p)
                    st.(key)=reshape(value,size(p));
                end
            end
        else
            T=T0*(p/p0).^((freestream.gamma-1)/freestream.gamma);
            st=struct('valid',isfinite(T)&T>=80&T<=2500,'temperature',T, ...
                'density',p./(freestream.R*T),'speed',sqrt(max(0,2*cp*(T0-T))));
        end
        atZero=ph==0;
        st.temperature(atZero)=T0; st.density(atZero)=rho0; st.speed(atZero)=0;
        st.pressure=p;
    end
end
function z=mu(T)
z=1.458e-6*T.^(3/2)./(T+110.4);
end
