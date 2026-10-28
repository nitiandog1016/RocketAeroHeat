function result = solveEntropyMatchedStagnation(freestream,geometry,wallTemperature,cfg)
%SOLVEENTROPYMATCHEDSTAGNATION 球锥鼻端熵层质量吞入与驻点层流剖面耦合。
% 模块作用：内部生成CR-112180物面关联激波表，求鼻端p0对流管质量
% 的导数，用可压缩驻点局部相似剖面近似闭合速度应变率与质量厚度。
% 这是鼻端局部极限；不求整个球锥、不接入生产外流或三维云图。
% 坐标为X主轴、鼻尖+X、零攻角速度-X；各量SI，Psi按轴对称每弧度。
%
% 推导：pe=ps-(ps-pinf)*(s/Rn)^2+O(s^4)，p0=ps+C*Psi+O(Psi^2)，
% Ue=a*s，Psi_BL=rhos*a*delta_m(a)*s^2。能量/熵恢复给出
%   a^2 = 2*(ps-pinf)/(rhos*Rn^2) + 2*C*a*delta_m(a)。
% C=4*A*(dp0/dcos^2Gamma)/(rhoinf*Uinf*Rn^2)，A取CR鼻端展开。
% 驻点local-similarity BVP(beta=.5)在Ue=0极限的无量纲剖面与a无关，
% 因而delta_m(a)正比a^(-1/2)。热流仍由实际耦合a下的同一剖面求得。
% 该式是已引用工程模型的局部推导，不是原文给出的通用热流关联式。
% 注意：质量吞入改变熵时，既定压力曲率与新a有有限的外缘动量残差。
% 本版BVP仍用A=-beta，只作显式容差内研究条件估计，不求解横向涡量。
% physicalBoundaryLayerClosed始终为false，available不等于源BL方程闭合。
% 条件：零攻角相切球锥、连续流假设、附着层流、薄层、同一热物性。
% 不含转捩/化学反应/稀薄流；未做几何全域或公开实验验证。
if nargin<4,cfg=struct;end
validateattributes(wallTemperature,{'numeric'},{'real','finite','scalar'});
opt=struct('maximumMassThicknessRatio',.1,'relativeTolerance',1e-9, ...
    'wallTemperatureSource','用户给定壁温边界；本模块不求壁温');
if isfield(cfg,'entropyStagnation')
    assert(isstruct(cfg.entropyStagnation)&&isscalar(cfg.entropyStagnation), ...
        'RocketAeroHeat:EntropyStagnation:Options','entropyStagnation必须是标量结构体。');
    names=fieldnames(cfg.entropyStagnation);
    for j=1:numel(names)
        assert(isfield(opt,names{j}),'RocketAeroHeat:EntropyStagnation:Options', ...
            '未知选项%s。',names{j});opt.(names{j})=cfg.entropyStagnation.(names{j});
    end
end
validateattributes(opt.maximumMassThicknessRatio,{'numeric'},{'scalar','finite','>',0,'<=',.2});
validateattributes(opt.relativeTolerance,{'numeric'},{'scalar','finite','>',0,'<=',1e-5});
assert((ischar(opt.wallTemperatureSource)||(isstring(opt.wallTemperatureSource)&& ...
    isscalar(opt.wallTemperatureSource)))&&~isempty(strtrim(char(opt.wallTemperatureSource))), ...
    'RocketAeroHeat:EntropyStagnation:Source','给定壁温必须有来源说明。');
result=struct('available',false,'status','鼻端质量/应变率耦合尚未完成', ...
    'method','CR-112180鼻端激波质量展开 + 同物性总压导数 + Levy-Lees驻点局部相似剖面', ...
    'scope','只求零攻角相切球锥驻点局部工程近似；保留外缘动量失配，未解熵层横向涡量；未验证全域或实验精度', ...
    'defaultSolverEnabled',false,'fullBodySolved',false,'experimentallyValidated',false, ...
    'physicalBoundaryLayerClosed',false,'entropyVorticityTermsSolved',false, ...
    'diagnosticOnly',true,'eligibleForSurfaceCoverage',false, ...
    'diagnosticAvailable',false,'conditionalEstimate',[], ...
    'localThermalMassClosureSolved',false,'continuumScreeningPassed',false, ...
    'wallTemperatureInput',wallTemperature,'wallTemperatureSource',char(opt.wallTemperatureSource), ...
    'velocityGradient',NaN,'fixedEntropyVelocityGradient',NaN,'massThickness',NaN, ...
    'heatFlux',NaN,'edgeTemperature',NaN,'edgePressure',NaN,'edgeDensity',NaN, ...
    'massClosureResidual',NaN,'derivativeRelativeUncertainty',NaN, ...
    'massThicknessRatio',NaN,'physicalThicknessRatio',NaN, ...
    'stagnationMomentumRelativeMismatch',NaN,'profile',[],'shockTable',[],'fieldMetadata',struct, ...
    'numericalOptions',opt);

%% 1. 内部压力、质量表及正激波驻点状态；任何无效状态不补值
table=solveAxisymmetricEntropyShockTable(freestream,geometry,cfg);
result.shockTable=table;
if ~table.available,result.status=['内部激波表不适用：' table.status];return;end
anchorCfg=cfg;anchorCfg.externalFlow.physicalPressureSource=table.physicalPressureSource;
anchorCfg.externalFlow.shockCrossingAngleSource='球鼻鼻端正激波穿越角Gamma=90度的解析极限';
anchor=recoverShockEntropyEdgeState(freestream,90,table.stagnationPressure,anchorCfg);
if ~anchor.available||anchor.speed~=0
    result.status='同热物性驻点状态未通过零速度/总焓核查';return;
end
ps=anchor.pressure;rs=anchor.density;Rn=geometry.noseRadius;
fixedA=sqrt(2*(ps-freestream.pressure)/(rs*Rn^2));
% 与当前生产链路使用同一动理论自由程/物理鼻半径筛选。它不是梯度Kn。
knLimit=.01;
if isfield(cfg,'heating')&&isfield(cfg.heating,'maxContinuumKnudsen')
    knLimit=cfg.heating.maxContinuumKnudsen;
end
validateattributes(knLimit,{'numeric'},{'scalar','finite','>',0,'<=',.01});
muInf=1.458e-6*freestream.temperature^1.5/(freestream.temperature+110.4);
muStag=1.458e-6*anchor.temperature^1.5/(anchor.temperature+110.4);
knInf=airMeanFreePath(freestream.pressure,freestream.temperature,muInf,freestream.R)/Rn;
knStag=airMeanFreePath(ps,anchor.temperature,muStag,freestream.R)/Rn;
result.continuumScreening=struct('freestreamKnudsen',knInf,'stagnationKnudsen',knStag, ...
    'maximumKnudsen',max(knInf,knStag),'limit',knLimit, ...
    'method','黏度动理论平均自由程 / 解析输入物理鼻端半径；不是梯度Kn');
if ~isfinite(max(knInf,knStag))||max(knInf,knStag)>=knLimit
    result.status='来流或驻点Kn超出连续流筛选；本模块不提供稀薄流解';return;
end
result.continuumScreeningPassed=true;

%% 2. 总压对穿越流管质量的鼻端导数；完全气体解析、NASA7三档收敛
if strcmp(table.thermalPropertyMode,'caloricallyPerfect')
    g=freestream.gamma;M=anchor.upstreamMachPhysical;Q=freestream.density*freestream.speed^2;
    p2=table.normalShockStaticPressure;D=2*Q/(g+1);
    B=ps*(D/(p2*(g-1))-2*g/((g-1)*((g-1)*M^2+2)));
    uncertainty=0;derivativeMethod='完全气体法向跃迁及完整总焓总压的解析鼻端导数';
else
    w=[1e-3 2.5e-4 6.25e-5];
    near=recoverShockEntropyEdgeState(freestream,acosd(sqrt(w)),freestream.pressure,cfg);
    if ~all(near.algebraicallyValid)
        result.status='同NASA7总压导数试探未通过热物性/守恒核查';return;
    end
    d=(near.totalPressure-ps)./w;
    previous=(4*d(2)-d(1))/3;B=(4*d(3)-d(2))/3;
    uncertainty=abs(B-previous)/abs(B);
    derivativeMethod='同NASA7完整总压对cos^2Gamma导数的两档Richardson外推';
end
if ~isfinite(B)||B<=0||uncertainty>1e-4
    result.status='鼻端总压导数不为正或未通过收敛核查';return;
end
C=4*table.noseExpansion.A*B/(freestream.density*freestream.speed*Rn^2);
result.totalPressureMassDerivative=C;result.totalPressureAngleDerivative=B;
result.derivativeRelativeUncertainty=uncertainty;result.derivativeMethod=derivativeMethod;

%% 3. 正激波固定熵应变率下先求剖面；只用于建立正确尺度的耦合方程
edge=struct('temperature',anchor.temperature,'pressure',ps,'density',rs,'speed',0, ...
    'totalEnthalpy',anchor.upstreamTotalEnthalpy,'sources',struct);
for key={'temperature','pressure','density','speed','totalEnthalpy'}
    edge.sources.(key{1})=['程序内部同热物性正激波驻点恢复；' table.physicalPressureSource];
end
station=struct('streamlineMetric',0,'xi',0,'beta',.5,'pressureGradient',0, ...
    'velocityGradient',fixedA,'stagnationVelocityGradient',fixedA,'sources',struct);
for key={'streamlineMetric','xi','beta','pressureGradient'}
    station.sources.(key{1})='零攻角球鼻解析驻点局部极限，h~s、U~a*s、xi~rho*mu*a*s^4/4';
end
station.sources.velocityGradient='同一内部压力曲率和正激波固定熵推得的初始应变率，非最终耦合值';
station.sources.stagnationVelocityGradient=station.sources.velocityGradient;
station.sources.wallTemperature=char(opt.wallTemperatureSource);
baseline=solveAxisymmetricLaminarMassProfile(freestream,edge,wallTemperature,station,cfg);
if ~baseline.available,result.status=['驻点局部剖面未解：' baseline.status];return;end

%% 4. 求真实正根；不以任意厚度或应变率替代质量吞入闭合
coefficient=2*C*baseline.massThickness*sqrt(fixedA);
equation=@(a)a.^2-fixedA^2-coefficient*sqrt(a);
upper=2*fixedA;expansions=0;
while equation(upper)<=0 && expansions<40,upper=2*upper;expansions=expansions+1;end
if ~isfinite(upper)||equation(upper)<=0
    result.status='正应变率质量闭合未找到可审查括区';return;
end
% fzero本身按当前x尺度应用TolX。不能再乘fixedA，否则高应变率下
% 等于重复放大容差，使重求剖面后的物理闭合残差无法通过。
matchedA=fzero(equation,[fixedA upper],optimset('TolX',opt.relativeTolerance));
station.velocityGradient=matchedA;station.stagnationVelocityGradient=matchedA;
station.sources.velocityGradient='内部激波总压/质量鼻端导数与同剖面质量厚度闭合求得的应变率';
station.sources.stagnationVelocityGradient=station.sources.velocityGradient;
profile=solveAxisymmetricLaminarMassProfile(freestream,edge,wallTemperature,station,cfg);
if ~profile.available,result.status=['耦合应变率下驻点剖面未解：' profile.status];return;end
residual=abs(matchedA^2-fixedA^2-2*C*matchedA*profile.massThickness)/matchedA^2;
ratio=profile.massThickness/Rn;
physicalRatio=profile.delta/Rn;
% s=0时两个梯度都为零，不能据此断言兼容。下式来自非零s的真正极限：
% A_lim=-(ps-pinf)/(rho_s*a^2*Rn^2)，beta=.5。
momentumMismatch=abs(matchedA^2-fixedA^2)/matchedA^2;
momentumLimit=profile.numericalOptions.maximumMomentumMismatch;
result.massClosureResidual=residual;result.massThicknessRatio=ratio;
result.physicalThicknessRatio=physicalRatio;
result.stagnationMomentumRelativeMismatch=momentumMismatch;
result.stagnationPressureParameter=-(ps-freestream.pressure)/(rs*matchedA^2*Rn^2);
result.stagnationVelocityParameter=.5;
result.momentumApproximation=struct('relativeMismatch',momentumMismatch,'limit',momentumLimit, ...
    'method','有限s压力/惯性梯度的鼻端相对失配；不是s=0的零除零检查', ...
    'scope','容差来自显式工程接受准则，非文献普适物理保证；熵层横向涡量项尚未求解');
if residual>10*opt.relativeTolerance||max(ratio,physicalRatio)>opt.maximumMassThicknessRatio
    result.status='鼻端质量残差或质量/几何薄层接受准则未通过';return;
end
% 即使动量不相容，也保留独立、明确带条件的研究诊断。主求解场保持
% NaN，不能为取得可用数值而放宽源剖面模块的默认兼容阈值。
result.diagnosticAvailable=true;result.localThermalMassClosureSolved=true;
result.conditionalEstimate=struct('velocityGradient',matchedA, ...
    'fixedEntropyVelocityGradient',fixedA,'massThickness',profile.massThickness, ...
    'heatFlux',profile.heatFlux,'edgeTemperature',anchor.temperature, ...
    'edgePressure',ps,'edgeDensity',rs,'profile',profile,'fixedEntropyProfile',baseline, ...
    'valueKind','conditionalEstimate','eligibleForSurfaceCoverage',false, ...
    'physicalClosureVerified',false,'source',result.method, ...
    'status','仅为条件代数诊断：未解熵层横向涡量，不能计入已求解面积');
if momentumMismatch>momentumLimit
    result.status='外缘动量鼻端极限失配，主求解场未解；条件代数诊断另存，不计入覆盖';return;
end
result.available=true;result.localThermalMassClosureSolved=true;
result.status='研究条件估计可用：代数质量/热力闭合，外缘动量失配已报告；不能计入表面求解覆盖';
result.velocityGradient=matchedA;result.fixedEntropyVelocityGradient=fixedA;
result.massThickness=profile.massThickness;result.heatFlux=profile.heatFlux;
result.edgeTemperature=anchor.temperature;result.edgePressure=ps;result.edgeDensity=rs;
result.profile=profile;result.fixedEntropyProfile=baseline;
result.sourceURLs={table.sourceURL,profile.sourceURL};
for key={'velocityGradient','massThickness','heatFlux','edgeTemperature','edgePressure','edgeDensity'}
    result.fieldMetadata.(key{1})=struct('valueKind','computed','valid',true, ...
        'source',result.method,'scope',result.scope);
end
for key={'velocityGradient','massThickness','heatFlux'}
    result.fieldMetadata.(key{1}).valueKind='conditionalEstimate';
    result.fieldMetadata.(key{1}).physicalClosureVerified=false;
    result.fieldMetadata.(key{1}).eligibleForSurfaceCoverage=false;
    result.fieldMetadata.(key{1}).status='研究条件估计；存在外缘动量失配，未闭合熵层横向涡量项';
end
result.fieldMetadata.wallTemperatureInput=struct('valueKind','prescribed','valid',true, ...
    'source',char(opt.wallTemperatureSource),'scope','给定壁温热边界，不是求得壁温');
end
