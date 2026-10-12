function edge = solveSphericalNoseEdge(mesh,freestream,cfg)
%SOLVESPHERICALNOSEEDGE 已验证球鼻的修正牛顿压力与驻点熵外缘近似。
% 模块作用：独立评估球鼻区域的 p、T、rho、U、M 和总状态，供后续
% 区域模型衔接研究。当前不接入默认外流链路，不覆盖已有网格或结果。
% 依据：NASA SP-8062 (1971), 4.2.1, 印刷页 30，钝鼻外缘可按
% 正激波后的熵膨胀到局部静压；该近似不能自动延伸到细长后体。
% 本模块以修正牛顿关系 p=p_inf+(p02-p_inf)*cos(theta)^2 给静压。
% 热力模型一致地计算 p02、T02 和局部状态；无高温定比热静默回退。
% 几何范围：重新核验球鼻、拟合范围内的完整面元、无特征/拓扑问题；
% 合成入射角 <=30 度，默认迎流极角 <=60 度（项目研究限制，非实验
% 验证界限）。不输出激波角/波面，不计算热流，不声称区域衔接已解。

%% 1. 输出契约：计算条件、算法来源和无解状态分别保存
if nargin<3, cfg=struct; end
if ~isfield(mesh,'localGeometry')
    error('RocketAeroHeat:SphereEdge:Geometry','需要当前网格的 localGeometry。');
end
nf=size(mesh.F,1); g=mesh.localGeometry;
validateattributes(g.center,{'numeric'},{'real','finite','size',[nf 3]});
validateattributes(g.normal,{'numeric'},{'real','finite','size',[nf 3]});
flow=double(freestream.unitFlowDirection(:)');
validateattributes(flow,{'numeric'},{'real','finite','size',[1 3]});
if abs(norm(flow)-1)>1e-8
    error('RocketAeroHeat:SphereEdge:Flow','来流方向必须为单位向量。');
end
limit=60;
if isfield(cfg,'externalFlow') && ...
        isfield(cfg.externalFlow,'sphericalNoseMaxPolarAngleDeg')
    limit=cfg.externalFlow.sphericalNoseMaxPolarAngleDeg;
end
validateattributes(limit,{'numeric'},{'real','finite','scalar','>',0,'<=',80});
[mode,thermalScope]=selectAirThermalModel(freestream,cfg);
blank=nan(nf,1);
edge=struct('available',false,'algebraicallyValid',false(nf,1), ...
    'geometryFace',false(nf,1),'scopeFace',false(nf,1), ...
    'pressure',blank,'temperature',blank,'density',blank,'speed',blank, ...
    'Mach',blank,'speedOfSound',blank,'totalPressure',blank, ...
    'totalTemperature',blank,'velocityVector',nan(nf,3), ...
    'incidenceNormal',nan(nf,3),'polarAngleDeg',blank, ...
    'energyResidual',blank,'entropyResidual',blank, ...
    'statusCode',ones(nf,1,'uint8'), ...
    'statusNames',{{'已计算（独立区域研究模型）','球鼻几何或工况不适用', ...
        '超出球鼻角度/几何范围','热物性或总状态不可用','切向流向退化'}}, ...
    'method',['球鼻修正牛顿压力 + 正激波驻点熵等熵恢复；' mode], ...
    'pressureMethod','解析球面径向入射角的修正牛顿压力', ...
    'stateMethod','同一热物性正激波总状态的等熵/总焓恢复', ...
    'normalSource','当前网格多尺度球面核验后的解析径向法向', ...
    'velocityMethod','来流在真实三角面切平面的投影，仅局部流向近似', ...
    'thermalPropertyMode',mode,'thermalScope',thermalScope, ...
    'sourceURL','https://ntrs.nasa.gov/api/citations/19710021703/downloads/19710021703.pdf', ...
    'sourceLocation','NASA SP-8062 (1971), 4.2.1, 印刷页 30 / PDF 页 36', ...
    'scope',sprintf(['独立研究模块；已核验球鼻，合成入射角<=30度，' ...
        '迎流极角<=%g度；仅拟合区域内完整光滑面元；不用于细长后体、' ...
        '分离、尾迹或波系干扰；未验证区域接口和实验精度'],limit), ...
    'maximumPolarAngleDeg',limit,'sphereCenter',[NaN NaN NaN], ...
    'noseRadius',NaN,'status','球鼻几何或工况尚未通过');

%% 2. 重新检查当前网格，不能用旧驻点缓存或给定半径代替几何证明
estimate=estimateNoseRadius(mesh);
edge.geometryEstimate=estimate;
if ~estimate.available
    edge.status=estimate.status; return;
end
incidence=acosd(max(-1,min(1,-flow(1))));
if freestream.Mach<=1 || incidence>30+1e-10
    edge.status='要求超声速且实际来流与 -X 轴夹角不超过 30 度'; return;
end
R=estimate.radius; V=double(mesh.V); F=double(mesh.F);
xTip=max(V(:,1)); tolerance=max(1e-10,1e-8*(xTip-min(V(:,1))));
tip=mean(V(abs(V(:,1)-xTip)<=tolerance,:),1);
origin=[xTip-R,tip(2:3)];
edge.sphereCenter=origin; edge.noseRadius=R;
radial=g.center-origin; radialLength=vecnorm(radial,2,2);
unitRadial=radial./radialLength;
mu=max(-1,min(1,unitRadial*(-flow')));
mu(abs(mu-1)<=8*eps)=1;
theta=acosd(mu);
vertexError=abs(vecnorm(V-origin,2,2)/R-1);
depth=xTip-V(:,1);
withinFit=depth<=estimate.outerDepthRadiusFraction*R+tolerance;
onFront=V(:,1)>=origin(1)-tolerance;
geometry=max(vertexError(F),[],2)<=0.03 & ...
    all(withinFit(F),2) & all(onFront(F),2) & ...
    sum(g.normal.*unitRadial,2)>=cosd(20) & ...
    ~g.featureFace & ~g.topologyIssue & g.isSidewall;
% 偏置驻点也必须位于已验证球面前缀；不得从拟合球心外推出缺失鼻端。
if 1-cosd(incidence)>estimate.outerDepthRadiusFraction+1e-8
    edge.status='偏置驻点超出已经验证的球面几何范围'; return;
end
edge.geometryFace=geometry;
edge.incidenceNormal(geometry,:)=unitRadial(geometry,:);
edge.polarAngleDeg(geometry)=theta(geometry);
inScope=geometry & theta<=limit+1e-10;
edge.scopeFace=inScope; edge.statusCode(:)=2;
edge.statusCode(inScope)=3;
if ~any(inScope)
    edge.status='没有满足球鼻几何及迎流极角范围的面元'; return;
end

%% 3. 正激波总状态与修正牛顿局部压力，保持同一热物性模型
if strcmp(mode,'thermallyPerfect')
    normal=solveThermallyPerfectNormalShock(freestream);
else
    normal=solveNormalShock(freestream);
end
edge.normalShock=normal;
if ~normal.available
    edge.status=normal.status; return;
end
ids=find(inScope); p0=normal.stagnationPressure;
T0=normal.stagnationTemperature;
p=freestream.pressure+(p0-freestream.pressure)*mu(ids).^2;
if strcmp(mode,'thermallyPerfect')
    state=recoverThermallyPerfectEdgeState(p,repmat(p0,size(p)),T0,freestream.R);
    valid=state.valid; T=state.temperature; rho=state.density;
    U=state.speed; M=state.Mach; sound=state.speedOfSound;
    [~,h0]=airThermalProperties(T0);
    [~,h]=airThermalProperties(T);
    energy=(h+0.5*U.^2-h0)/max(1,abs(h0));
    entropy=state.entropyResidual;
else
    cp=freestream.gamma*freestream.R/(freestream.gamma-1);
    T=T0*(p/p0).^((freestream.gamma-1)/freestream.gamma);
    U=sqrt(max(0,2*cp*(T0-T)));
    rho=p./(freestream.R*T);
    sound=sqrt(freestream.gamma*freestream.R*T); M=U./sound;
    energy=(cp*T+0.5*U.^2-cp*T0)/(cp*T0);
    entropy=cp*log(T0./T)-freestream.R*log(p0./p);
    valid=isfinite(T) & T>0 & isfinite(U) & p>0 & p<=p0;
end

%% 4. 实际三角面切向速度及逐面元状态；无解字段保持 NaN
% p=p0 的解析驻点速度严格为零，不能把熵反算二分舍入误差开平方
% 成非零速度。零速度是求得的极限值，与域外面元的 NaN 明确区分。
atStagnation=valid & p==p0;
T(atStagnation)=T0; U(atStagnation)=0; M(atStagnation)=0;
rho(atStagnation)=p0/(freestream.R*T0);
energy(atStagnation)=0; entropy(atStagnation)=0;
tangent=flow-sum(g.normal(ids,:).*flow,2).*g.normal(ids,:);
tangentLength=vecnorm(tangent,2,2);
degenerate=valid & tangentLength<=1e-10 & U>1e-8;
edge.statusCode(ids(degenerate))=4;
valid=valid & ~degenerate;
keep=ids(valid);
edge.pressure(keep)=p(valid); edge.temperature(keep)=T(valid);
edge.density(keep)=rho(valid); edge.speed(keep)=U(valid);
edge.Mach(keep)=M(valid); edge.speedOfSound(keep)=sound(valid);
edge.totalPressure(keep)=p0; edge.totalTemperature(keep)=T0;
edge.energyResidual(keep)=energy(valid); edge.entropyResidual(keep)=entropy(valid);
velocity=zeros(numel(ids),3);
moving=valid & tangentLength>1e-10;
velocity(moving,:)=U(moving).*tangent(moving,:)./tangentLength(moving);
edge.velocityVector(keep,:)=velocity(valid,:);
edge.algebraicallyValid(keep)=true; edge.statusCode(keep)=0;
edge.available=any(edge.algebraicallyValid);
edge.status=sprintf('独立球鼻外缘模型计算 %d 面元；尚未接入默认链路或验证区域接口',numel(keep));
end
