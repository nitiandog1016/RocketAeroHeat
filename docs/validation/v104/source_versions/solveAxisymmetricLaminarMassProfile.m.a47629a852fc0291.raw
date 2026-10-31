function result = solveAxisymmetricLaminarMassProfile(freestream,edge,wallTemperature,station,cfg)
%SOLVEAXISYMMETRICLAMINARMASSPROFILE 独立轴对称层流局部相似剖面。
% 模块作用：解同一速度/总焓剖面，直接积分物理动量、位移及质量
% 厚度，并求壁面热流、剪应力。不是 Thwaites 厚度换算，不采用5.27。
% 来源：NASA TP-3271(1992)，式6–7、17、25、54、56–61；本版把
% 沿 xi 的非相似项置零，属于 local similarity 工程近似。
% 不预测转捩、分离/再附、横流、涡量相互作用或实验精度。
% 暖壁有利压梯下的附着支允许速度略高于外缘再返回外缘；
% 速度全域单调不是附着条件。来源NACA TN-3325印页17–18、22。
%
% edge 标量字段：temperature/pressure/density/speed/totalEnthalpy，
% edge.sources 同名来源文字。station 字段：streamlineMetric [m]、
% xi [kg^2/s^2]、beta、pressureGradient [Pa/m]、velocityGradient [1/s]；
% station.sources 保存同名来源及 wallTemperature 来源。普通站 U,h,xi>0。
% 解析驻点允许 U=h=xi=0，另给 stagnationVelocityGradient=a>0；beta=.5，
% K=sqrt(nu_e/(2a))，不以任意 epsilon 启动。Twall [K] 为给定边界。
% cfg.massProfile 可设置求解/外缘/动量兼容容差；cfg.heating 中的
% airThermalModel 沿用项目明确的 NASA7 或 caloricallyPerfect 选择。

%% 1. 校验字段和来源；物理域外保留未解，不用输入值冒充求解结果
if nargin<5, cfg=struct; end
checkStruct(freestream,'freestream'); checkStruct(edge,'edge');
checkStruct(station,'station'); checkStruct(cfg,'cfg');
requireNumbers(freestream,{'R','gamma'},'freestream');
requireNumbers(edge,{'temperature','pressure','density','speed','totalEnthalpy'},'edge');
requireNumbers(station,{'streamlineMetric','xi','beta','pressureGradient','velocityGradient'},'station');
validateattributes(wallTemperature,{'numeric'},{'real','finite','scalar'});
requireSources(edge,{'temperature','pressure','density','speed','totalEnthalpy'},'edge');
requireSources(station,{'streamlineMetric','xi','beta','pressureGradient', ...
    'velocityGradient','wallTemperature'},'station');
[mode,thermalScope]=selectAirThermalModel(freestream,cfg);
thermalScope=standaloneAirThermalScope(mode);
opt=options(cfg);
result=blankResult(edge,wallTemperature,station,mode,thermalScope,opt);
R=freestream.R; cpCal=freestream.gamma*R/(freestream.gamma-1);
Te=edge.temperature; pe=edge.pressure; rhoe=edge.density; Ue=edge.speed;
metric=station.streamlineMetric; xi=station.xi; beta=station.beta;
if R<=0 || freestream.gamma<=1 || pe<=0 || rhoe<=0 || Ue<0 || ...
        Te<80 || Te>2500 || wallTemperature<80 || wallTemperature>2500 || ...
        beta<0 || beta>1 || metric<0 || xi<0
    result=refuse(result,1,'输入超出本版非负 beta、附着层流及80~2500 K研究范围'); return;
end
if abs(rhoe*R*Te/pe-1)>1e-6
    result=refuse(result,1,'外缘密度与静压、静温、气体常数不一致'); return;
end

%% 2. 统一焓零点及输运性质，核对总焓与沿程压力/速度梯度
% NASA7 的原焓零点任意。数值 g 使用正的移位总焓；输出还原原焓。
% 此平移只改变数值变量，不改变热流、能量或物理温度。
if strcmp(mode,'thermallyPerfect')
    [cp80,h80]=airThermalProperties(80);
    offset=h80-cp80*80;
    [~,he]=airThermalProperties(Te); [~,hw]=airThermalProperties(wallTemperature);
else
    offset=0; he=cpCal*Te; hw=cpCal*wallTemperature;
end
Hphysical=edge.totalEnthalpy; H0=Hphysical-offset;
energyMismatch=abs(Hphysical-he-Ue^2/2)/max([abs(he-offset),Ue^2/2,1]);
result.inputEnergyResidual=energyMismatch; result.enthalpyOffset=offset;
if H0<=0 || energyMismatch>opt.inputEnergyTolerance
    result=refuse(result,1,'外缘总焓与同物性静温、速度不一致'); return;
end
muEdge=sutherland(Te); nuEdge=muEdge/rhoe;
stagnation=(Ue==0 && metric==0 && xi==0);
if stagnation
    requireNumbers(station,{'stagnationVelocityGradient'},'station');
    requireSources(station,{'stagnationVelocityGradient'},'station');
    a=station.stagnationVelocityGradient;
    if a<=0 || abs(beta-.5)>1e-10 || abs(station.velocityGradient/a-1)>1e-6 || ...
            station.pressureGradient~=0
        result=refuse(result,1,'解析轴对称驻点需要a>0、beta=.5及一致的驻点梯度'); return;
    end
    K=sqrt(nuEdge/(2*a)); expectedBeta=.5;
    mismatch=0; rawResidual=0; scale=0;
    result.assumedStagnationPressureCurvature=-rhoe*a^2;
    result.stagnationMomentumScope=['仅验证驻点点值U=pe_s=0；beta=.5的BVP另假设' ...
        'A=-.5，即近驻点pe_ss=-rhoe*a^2；本接口未独立核对该二阶曲率。' ...
        '不得由点值零残差声称熵层吞入后的近驻点动量兼容'];
elseif Ue>0 && metric>0 && xi>0
    K=sqrt(2*xi)/(rhoe*Ue*metric);
    expectedBeta=2*xi*station.velocityGradient/(rhoe*muEdge*Ue^2*metric^2);
    velocityTerm=rhoe*Ue*station.velocityGradient;
    rawResidual=station.pressureGradient+velocityTerm;
    scale=max(abs(station.pressureGradient),abs(velocityTerm));
    if scale==0, mismatch=0; else, mismatch=abs(rawResidual)/scale; end
else
    result=refuse(result,1,'普通站U、h、xi均须为正；只接受完整解析驻点零极限'); return;
end
result.stagnationLimit=stagnation;
result.thicknessScale=K; result.expectedBeta=expectedBeta;
result.inputMomentumResidual=rawResidual;
result.inputMomentumRelativeResidual=mismatch;
result.inputMomentumDenominator=scale;
result.inputMomentumDenominatorDefinition='max(abs(dpe/ds),abs(rhoe*Ue*dUe/ds))；两项均零时相对残差为零';
if abs(beta-expectedBeta)>opt.gradientConsistencyTolerance*max(1,abs(expectedBeta))
    result=refuse(result,7,'beta 与同一 xi、度量及外缘速度梯度不一致'); return;
end
if mismatch>opt.maximumMomentumMismatch
    result=refuse(result,2,'压力与速度梯度不满足本版动量兼容容差'); return;
end
result.momentumClosure=['局部相似方程采用A=-beta；忽略外缘熵层横向涡量贡献；' ...
    '容许残差阈值是显式工程接受准则，不是文献普适范围'];
gw=(hw-offset)/H0; Pr=opt.Prandtl;
trialExtension=false;

%% 3. 五状态BVP：f、F、m=l F_eta、g、n=(l/Pr)g_eta+omega
% NASA7反焓只在数值迭代试探阶段延拓；最终全剖面严格检查原物性域。
% 拒绝超域终解，绝不把试探延拓温度当作最终求解值。
etaMax=opt.etaMax;
try
    initial=bvpinit(linspace(0,etaMax,81),@guess);
    if isfield(opt,'initialProfile') && ~isempty(opt.initialProfile)
        ip=opt.initialProfile;
        if ~isstruct(ip) || ~isfield(ip,'eta') || ~isfield(ip,'states') || ...
                size(ip.states,1)~=5 || numel(ip.eta)~=size(ip.states,2) || ...
                any(~isfinite(ip.states(:))) || any(diff(ip.eta)<=0)
            error('RocketAeroHeat:LaminarMassProfile:Initial','初始剖面须为有限、递增eta及5行states。');
        end
        initial=bvpinit(linspace(0,etaMax,81), ...
            @(eta)interp1(ip.eta(:),ip.states',eta,'pchip','extrap')');
    end
    bopt=bvpset('RelTol',opt.relativeTolerance,'AbsTol',opt.absoluteTolerance, ...
        'NMax',opt.maximumMeshPoints,'Stats','off');
    sol=bvp4c(@ode,@boundary,initial,bopt);
catch exception
    if strcmp(exception.identifier,'RocketAeroHeat:LaminarMassProfile:Initial'), rethrow(exception); end
    result.numericalMessage=exception.message;
    result=refuse(result,4,'局部相似剖面BVP未取得可审查数值解'); return;
end
result.trialPropertyExtensionUsed=trialExtension;

%% 4. 终解逐点核查及独立外缘判定；etaMax仅为数值截断
eta=linspace(0,etaMax,opt.auditPoints);
[Y,Yp]=deval(sol,eta);
[T,l,D,physical]=properties(Y);
if ~all(physical) || any(~isfinite(Y(:)))
    result=refuse(result,3,'最终剖面超出原热物性范围，数值试探延拓不计求解'); return;
end
rhs=ode(eta,Y);
odeResidual=max(abs(Yp-rhs)./(1+abs(rhs)),[],'all');
boundaryResidual=max(abs(boundary(Y(:,1),Y(:,end))));
F=Y(2,:); m=Y(3,:); g=Y(4,:); n=Y(5,:);
omega=(1-1/Pr)*(Ue^2/H0)*F.*m;
dF=m./l; dg=Pr*(n-omega)./l;
result.odeResidual=odeResidual; result.boundaryResidual=boundaryResidual;
% 模块作用：记录真实速度支形状，区分无逆流附着支与数值拒绝。
% 有利压梯作用于温度较高、密度较小的近外缘流层，可使u/Ue>1；
% TP-3271式19及原边界允许随后从上方返回外缘。保留这些真实值，
% 不裁剪速度、不把负的近外缘导数等同于负壁面剪应力或分离。
result.velocityAudit=struct('minimumRatio',min(F),'maximumRatio',max(F), ...
    'overshootRatio',max(0,max(F)-1),'minimumEtaDerivative',min(dF), ...
    'wallMomentumSlope',m(1),'noBackflow',min(F)>=-opt.monotonicityTolerance, ...
    'positiveWallShearSlope',m(1)>0, ...
    'auditDomain','完整数值剖面0<=eta<=etaMax；边缘另由多场量独立核查', ...
    'scope','原始支形状诊断；过冲大小本身不是热流实验精度或附着物理充分证明');
if odeResidual>opt.profileResidualTolerance || boundaryResidual>1e-7 || ...
        m(1)<=0 || min(F)<-opt.monotonicityTolerance
    result=refuse(result,6,'剖面残差、正壁面剪应力或无逆流附着支检查失败'); return;
end
vDefect=abs(1-F)/opt.edgeVelocityTolerance;
hDefect=abs(1-g)/opt.edgeEnthalpyTolerance;
% 总焓受焓零点和大动能影响，不能单独证明静温/密度已回到外缘。
% 同时约束参考零点无关的T/Te与rho/rhoe，避免高焓工况提前截断。
tDefect=abs(T/Te-1)/opt.edgeTemperatureTolerance;
rhoDefect=abs(Te./T-1)/opt.edgeDensityTolerance;
shearDefect=abs(dF)/max(max(abs(dF)),eps)/opt.edgeGradientTolerance;
% 零热流极限中g_eta可小于BVP绝对精度；只作相对归一化会放大噪声，
% 令物理相同的驻点与近驻点得到不同质量边缘。采用显式混合数值容差；
% 绝对阈值不是驻点启动epsilon，也不改变速度/温度/密度边缘要求。
heatGradientMaximum=max(abs(dg));
heatGradientThreshold=max(opt.edgeGradientTolerance*heatGradientMaximum, ...
    opt.edgeAbsoluteGradientTolerance);
heatDefect=abs(dg)/heatGradientThreshold;
defect=max([vDefect;hDefect;tDefect;rhoDefect;shearDefect;heatDefect],[],1);
lastBad=find(defect>1,1,'last');
if isempty(lastBad) || lastBad>=numel(eta) || eta(lastBad+1)>etaMax-opt.edgeBuffer
    result=refuse(result,5,'数值外界距真实有限边缘不足，增大etaMax并检查收敛'); return;
end
etaEdge=fzero(@edgeDefect,[eta(lastBad),eta(lastBad+1)]);
eY=deval(sol,etaEdge); [eT,el,~,~]=properties(eY);
em=eY(3); eg=eY(4); en=eY(5); eF=eY(2);
eOmega=(1-1/Pr)*(Ue^2/H0)*eF*em;
absoluteHeatGradient=abs(Pr*(en-eOmega)/el);
edgeChecks=struct('velocityDefect',abs(1-eF),'totalEnthalpyDefect',abs(1-eg), ...
    'relativeTemperatureDefect',abs(eT/Te-1),'relativeDensityDefect',abs(Te/eT-1), ...
    'relativeShearGradient',abs(em/el)/max(max(abs(dF)),eps), ...
    'relativeHeatGradient',0,'absoluteHeatGradient',absoluteHeatGradient, ...
    'heatGradientMaximum',heatGradientMaximum,'heatGradientThreshold',heatGradientThreshold, ...
    'heatGradientCriterionDefect',absoluteHeatGradient/heatGradientThreshold, ...
    'heatGradientCriterion','abs(g_eta)<=max(relativeTol*max(abs(g_eta)),absoluteTol)', ...
    'temperature',eT);
if heatGradientMaximum>0
    edgeChecks.relativeHeatGradient=absoluteHeatGradient/heatGradientMaximum;
end

%% 5. 用同一有限剖面独立积分厚度、壁面热流及剪应力
% 质量积分的上限明确为etaEdge；不把无限域、etaMax或旧theta当质量边缘。
nodes=linspace(0,etaEdge,opt.integrationPoints);
states=deval(sol,nodes); [temp,ell,ratio,~]=properties(states);
massIntegral=integral(@massIntegrand,0,etaEdge,'RelTol',1e-8,'AbsTol',1e-10);
momentumIntegral=integral(@momentumIntegrand,0,etaEdge,'RelTol',1e-8,'AbsTol',1e-10);
deltaIntegral=integral(@thicknessIntegrand,0,etaEdge,'RelTol',1e-8,'AbsTol',1e-10);
displacementIntegral=deltaIntegral-massIntegral;
delta=K*deltaIntegral; ds=K*displacementIntegral;
theta=K*momentumIntegral; mass=K*massIntegral;
lwall=l(1); qw=lwall*muEdge*H0/(Pr*K)*dg(1);
tau=lwall*muEdge*Ue/K*dF(1);
if any(~isfinite([delta ds theta mass qw tau])) || delta<=0 || theta<=0 || mass<=0
    result=refuse(result,6,'有限剖面厚度或壁面场量的物理积分检查失败'); return;
end
result.available=true; result.algebraicallyValid=true; result.statusCode=uint8(0);
result.status='已计算：独立局部相似层流近似；未验证几何耦合或实验精度';
result.delta=delta; result.deltaStar=ds; result.thetaPhysical=theta;
result.massThickness=mass; result.heatFlux=qw; result.wallShearStress=tau;
result.massFluxPerMetric=rhoe*Ue*metric*mass;
result.wallTemperatureInput=wallTemperature;
result.etaEdge=etaEdge; result.edgeChecks=edgeChecks;
result.profile=struct('eta',nodes,'states',states,'temperature',temp, ...
    'speed',Ue*states(2,:),'density',rhoe./ratio, ...
    'totalEnthalpy',offset+H0*states(4,:), ...
    'viscosity',sutherland(temp),'physicalY',K*cumtrapz(nodes,ratio), ...
    'velocityDerivative',states(3,:)./ell);
result.numericalProfile=struct('eta',sol.x,'states',sol.y);
result.dimensionlessIntegrals=struct('mass',massIntegral,'momentum',momentumIntegral, ...
    'physicalThickness',deltaIntegral,'displacement',displacementIntegral);
result.wallDimensionlessVelocitySlope=dF(1);
result.wallDimensionlessEnthalpySlope=dg(1);
% BVP的AbsTol/RelTol是数值尺度而非严格误差界；极小热流不得宣称
% 具有相对精度。只标明壁面g_eta是否显著高于其状态容差换算尺度。
result.wallHeatGradientNumericalScale=Pr/lwall* ...
    (opt.absoluteTolerance+opt.relativeTolerance*abs(n(1)));
result.heatFluxNumericallyResolved=abs(dg(1))> ...
    opt.heatGradientResolutionFactor*result.wallHeatGradientNumericalScale;
result.heatFluxResolutionScope=['依据BVP状态容差尺度的数值分辨标记，非严格误差界或实验精度；' ...
    'false时不得声称极小热流具可靠相对精度，原热流与导数不裁剪'];
result.massDefinitionResidual=abs(delta-ds-mass)/max(delta,mass);
result.fieldMetadata=syncMetadata(result.fieldMetadata,true);
result.fieldMetadata.heatFlux.numericallyResolved=result.heatFluxNumericallyResolved;
result.fieldMetadata.heatFlux.resolutionScope=result.heatFluxResolutionScope;

    function dy=ode(~,y)
        [~,lc,dc,~]=properties(y);
        ff=y(1,:); vv=y(2,:); mm=y(3,:); nn=y(5,:);
        ww=(1-1/Pr)*(Ue^2/H0)*vv.*mm;
        vp=mm./lc; gp=Pr*(nn-ww)./lc;
        dy=[vv;vp;-ff.*vp+beta*(vv.^2-dc);gp;-ff.*gp];
    end
    function residual=boundary(a,b)
        residual=[a(1);a(2);a(4)-gw;b(2)-1;b(4)-1];
    end
    function y=guess(etaGuess)
        v=tanh(etaGuess); vp=1-v^2;
        gg=gw+(1-gw)*v;
        raw=[log(cosh(etaGuess));v;vp;gg;0];
        [~,lc,~,~]=properties(raw);
        mm=lc*vp;
        ww=(1-1/Pr)*(Ue^2/H0)*v*mm;
        y=[raw(1);v;mm;gg;(lc/Pr)*(1-gw)*vp+ww];
    end
    function [t,lc,dc,valid]=properties(y)
        ht=H0*y(4,:)-Ue^2*y(2,:).^2/2+offset;
        if strcmp(mode,'caloricallyPerfect')
            rawT=ht/cpCal;
        else
            rawT=inverseNASAEnthalpy(ht);
        end
        valid=isfinite(rawT) & rawT>=80-1e-7 & rawT<=2500+1e-7;
        if any(~valid), trialExtension=true; end
        % 只为Newton试探构造正输运系数；终解以rawT原值严格拒绝域外。
        t=max(1,rawT);
        dc=t/Te;
        lc=(Te./t).*sutherland(t)/muEdge;
    end
    function d=edgeDefect(et)
        yy=deval(sol,et); [tt,lc,~,~]=properties(yy);
        ww=(1-1/Pr)*(Ue^2/H0)*yy(2)*yy(3);
        hh=abs(Pr*(yy(5)-ww)/lc)/heatGradientThreshold;
        d=max([abs(1-yy(2))/opt.edgeVelocityTolerance, ...
            abs(1-yy(4))/opt.edgeEnthalpyTolerance, ...
            abs(tt/Te-1)/opt.edgeTemperatureTolerance, ...
            abs(Te/tt-1)/opt.edgeDensityTolerance, ...
            abs(yy(3)/lc)/max(max(abs(dF)),eps)/opt.edgeGradientTolerance,hh])-1;
    end
    function v=massIntegrand(et)
        yy=deval(sol,et); v=yy(2,:);
    end
    function v=momentumIntegrand(et)
        yy=deval(sol,et); v=yy(2,:).*(1-yy(2,:));
    end
    function v=thicknessIntegrand(et)
        yy=deval(sol,et); [~,~,v,~]=properties(yy);
    end
end

function requireNumbers(s,names,label)
% 模块作用：结构错误直接报告，数值物理域由调用处保存拒绝状态。
for k=1:numel(names)
    key=names{k};
    if ~isfield(s,key), error('RocketAeroHeat:LaminarMassProfile:Input','缺少%s.%s。',label,key); end
    validateattributes(s.(key),{'numeric'},{'real','finite','scalar'},mfilename,[label '.' key]);
end
end
function checkStruct(s,label)
if ~isstruct(s) || ~isscalar(s)
    error('RocketAeroHeat:LaminarMassProfile:Input','%s必须为标量结构体。',label);
end
end
function requireSources(s,names,label)
if ~isfield(s,'sources') || ~isstruct(s.sources) || ~isscalar(s.sources)
    error('RocketAeroHeat:LaminarMassProfile:Source','%s缺少sources逐字段来源。',label);
end
for k=1:numel(names)
    key=names{k};
    if ~isfield(s.sources,key) || ...
            ~(ischar(s.sources.(key)) || (isstring(s.sources.(key)) && isscalar(s.sources.(key)))) || ...
            isempty(strtrim(char(s.sources.(key))))
        error('RocketAeroHeat:LaminarMassProfile:Source','缺少%s.sources.%s的有效来源。',label,key);
    end
end
end
function o=options(cfg)
% 模块作用：所有人为数值/适用阈值显式保存，便于研究收敛与复核。
% monotonicityTolerance保留旧配置名称以兼容历史输入；现只用于
% 最小速度的数值逆流容差，不再限制F<=1或全域F_eta>=0。
o=struct('Prandtl',.72,'etaMax',12,'relativeTolerance',1e-7, ...
    'absoluteTolerance',1e-9,'maximumMeshPoints',12000,'auditPoints',1201, ...
    'integrationPoints',601,'edgeVelocityTolerance',.005, ...
    'edgeEnthalpyTolerance',.005,'edgeGradientTolerance',.005, ...
    'edgeAbsoluteGradientTolerance',1e-8, ...
    'edgeTemperatureTolerance',.005,'edgeDensityTolerance',.005, ...
    'edgeBuffer',1,'maximumMomentumMismatch',1e-6, ...
    'gradientConsistencyTolerance',1e-6,'inputEnergyTolerance',1e-6, ...
    'profileResidualTolerance',2e-5,'monotonicityTolerance',1e-6);
o.heatGradientResolutionFactor=10;
if isfield(cfg,'heating') && isfield(cfg.heating,'Prandtl'), o.Prandtl=cfg.heating.Prandtl; end
if isfield(cfg,'massProfile')
    checkStruct(cfg.massProfile,'cfg.massProfile');
    keys=fieldnames(cfg.massProfile);
    for k=1:numel(keys)
        key=keys{k};
        if ~isfield(o,key) && ~strcmp(key,'initialProfile')
            error('RocketAeroHeat:LaminarMassProfile:Option','未知massProfile选项%s。',key);
        end
        o.(key)=cfg.massProfile.(key);
    end
end
keys=fieldnames(o);
for k=1:numel(keys)
    if strcmp(keys{k},'initialProfile'), continue; end
    validateattributes(o.(keys{k}),{'numeric'},{'real','finite','scalar','positive'});
end
if o.Prandtl>2 || o.etaMax<4 || o.etaMax>60 || o.edgeBuffer>=o.etaMax || ...
        o.maximumMomentumMismatch>1 || o.edgeVelocityTolerance>=.1 || ...
        o.edgeEnthalpyTolerance>=.1 || o.edgeGradientTolerance>=.1 || ...
        o.edgeTemperatureTolerance>=.1 || o.edgeDensityTolerance>=.1 || ...
        o.edgeAbsoluteGradientTolerance>=.1
    error('RocketAeroHeat:LaminarMassProfile:Option','数值选项或外缘/兼容容差超出研究模块范围。');
end
for key={'maximumMeshPoints','auditPoints','integrationPoints'}
    validateattributes(o.(key{1}),{'numeric'},{'integer','>=',100});
end
end
function r=blankResult(edge,Tw,station,mode,thermalScope,opt)
% 模块作用：输入/输出、来源与是否真正求得分开；无解字段保持NaN。
r=struct('available',false,'algebraicallyValid',false,'statusCode',uint8(1), ...
    'status','尚无通过核查的局部相似解', ...
    'statusNames',{{'已计算的局部相似近似','输入范围或同物性状态不一致', ...
        '外缘动量梯度不兼容','终解热物性域外','剖面BVP失败', ...
        '有限外缘未解析','附着剖面或残差检查失败','beta与输入梯度不一致'}}, ...
    'method','轴对称Levy-Lees局部相似层流速度/总焓BVP + 同剖面物理积分', ...
    'scope',['独立local similarity；0<=beta<=1、正壁面剪应力及无逆流附着支；' ...
        '暖壁有利压梯支允许返回外缘的速度过冲；' ...
        '无转捩、分离再附或横向涡量交互；未验证真实几何耦合或实验精度'], ...
    'sourceURL','https://ntrs.nasa.gov/api/citations/19930006245/downloads/19930006245.pdf', ...
    'sourceEquations','式6–7、17、25、54、56–61；沿xi非相似项在本版置零', ...
    'localSimilarity',true,'nonSimilarTermsSolved',false,'geometryVerified',false, ...
    'massMatchingSolved',false,'thermalPropertyMode',mode,'thermalScope',thermalScope, ...
    'stagnationPressureCurvatureVerified',false, ...
    'edgeInput',edge,'stationInput',station,'wallTemperatureInput',Tw, ...
    'numericalOptions',opt,'profile',[],'numericalProfile',[], ...
    'fieldMetadata',struct,'trialPropertyExtensionUsed',false);
r.velocityAudit=[];
r.attachedBranchCriterion=['原物性域、ODE/边界残差、正壁面剪应力、无逆流及独立有限外缘；' ...
    '不要求速度全域单调，不裁剪热扩散/有利压梯引起的过冲；非实验准确性声明'];
r.velocityOvershootSourceURL='https://ntrs.nasa.gov/api/citations/19930084102/downloads/19930084102.pdf';
r.velocityOvershootSourceScope=['NACA TN-3325印页17–18及22说明暖壁有利压梯的粘性支可过冲；' ...
    '该报告使用Pr=1、线性黏性温度关系，不作为本NASA7/Pr0.72热流定量基准'];
r.heatFluxNumericallyResolved=false;
fields={'delta','deltaStar','thetaPhysical','massThickness','heatFlux','wallShearStress', ...
    'massFluxPerMetric','etaEdge','thicknessScale','inputEnergyResidual', ...
    'inputMomentumResidual','inputMomentumRelativeResidual','inputMomentumDenominator', ...
    'odeResidual','boundaryResidual','wallHeatGradientNumericalScale', ...
    'assumedStagnationPressureCurvature'};
for k=1:numel(fields), r.(fields{k})=NaN; end
for key={'delta','deltaStar','thetaPhysical','massThickness','heatFlux','wallShearStress','massFluxPerMetric'}
    r.fieldMetadata.(key{1})=struct('valueKind','computed','valid',false, ...
        'source',r.method,'scope',r.scope);
end
r.fieldMetadata.wallTemperatureInput=struct('valueKind','prescribed','valid',true, ...
    'source',char(station.sources.wallTemperature),'scope','给定热边界；本模块不求壁温');
for key={'pressure','temperature','density','speed','totalEnthalpy'}
    r.fieldMetadata.(key{1})=struct('valueKind','input','valid',true, ...
        'source',char(edge.sources.(key{1})),'scope','输入外缘状态，未由本剖面模块求得');
end
end
function r=refuse(r,code,message)
r.statusCode=uint8(code); r.status=message;
r.fieldMetadata=syncMetadata(r.fieldMetadata,false);
end
function m=syncMetadata(m,valid)
keys=fieldnames(m);
for k=1:numel(keys)
    if strcmp(m.(keys{k}).valueKind,'computed'), m.(keys{k}).valid=valid; end
end
end
function mu=sutherland(T)
mu=1.458e-6*T.^(3/2)./(T+110.4);
end
function T=inverseNASAEnthalpy(h)
% 模块作用：缓存原NASA7焓表并反解；最终用原函数Newton细化。
% 表外仅供BVP试探做线性延拓，不产生可用输出状态。
persistent hGrid tGrid cpLo cpHi hLo hHi
if isempty(hGrid)
    tGrid=(80:1:2500)';
    [~,hGrid]=airThermalProperties(tGrid);
    [cpLo,hLo]=airThermalProperties(80); [cpHi,hHi]=airThermalProperties(2500);
end
T=interp1(hGrid,tGrid,h,'pchip');
low=h<hLo; high=h>hHi;
T(low)=80+(h(low)-hLo)/cpLo;
T(high)=2500+(h(high)-hHi)/cpHi;
inside=~low & ~high & isfinite(h);
if any(inside)
    t=T(inside);
    for k=1:3
        [cp,hh]=airThermalProperties(t); t=t+(h(inside)-hh)./cp;
    end
    T(inside)=t;
end
end
