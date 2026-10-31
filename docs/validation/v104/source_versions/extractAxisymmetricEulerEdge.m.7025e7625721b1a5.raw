function edge = extractAxisymmetricEulerEdge(result,options)
%EXTRACTAXISYMMETRICEULEREDGE 从内部 Euler 场提取近壁外缘候选与误差诊断。
% 模块用途：为后续边界层接口准备沿母线的 p、T、rho、U 和 dU/ds。
% 一阶解使用首层单元近似；线性重构解使用 r 加权面位置的内侧迹，
% 然后把速度投影到几何切向，保留被投影掉的法向分量作为诊断。
% 这不是已验证的边界层外缘解；不借用自由来流总压、不强制驻点速度，
% 不混入反射壁面通量压力反算温度。原始单元及壁压差异分别保存。
% 坐标 s 从 +X 鼻端沿母线向尾部增大；输入和输出均为 SI 单位。
if nargin<2,options=struct;end
cornerAngle=5;
if isfield(options,'cornerAngleDeg'),cornerAngle=options.cornerAngleDeg;end
validateattributes(cornerAngle,{'numeric'},{'real','finite','scalar','>',0,'<',90});
g=result.grid;fs=result.freestream;
if ~strcmp(result.thermodynamicModel,'caloricallyPerfect')
    error('RocketAeroHeat:EulerEdge:Thermodynamics','当前提取只支持定比热 Euler 场。');
end
for name={'density','pressure','velocityX','velocityR'}
    validateattributes(result.(name{1}),{'numeric'},{'real','finite','size',[g.nx,g.nr]});
end
if any(result.density<=0,'all') || any(result.pressure<=0,'all')
    error('RocketAeroHeat:EulerEdge:State','Euler 场密度和压力必须为正。');
end
validateattributes(fs.gamma,{'numeric'},{'real','finite','scalar','>',1});
validateattributes(fs.R,{'numeric'},{'real','finite','scalar','positive'});
for name={'alpha','beta'}
    if isfield(fs,name{1}) && (~isfinite(fs.(name{1})) || abs(fs.(name{1}))>1e-10)
        error('RocketAeroHeat:EulerEdge:Incidence','轴对称外缘提取要求零攻角、零侧滑。');
    end
end

%% 1. 从真实网格壁面计算弧长和下游切向；鼻前轴线不作为机体表面
dx=diff(g.X(:,1));dr=diff(g.R(:,1));ds=hypot(dx,dr);
wall=g.bodyWall & g.wallR>0;
index=flipud(find(wall));n=numel(index);
if n<3 || any(diff(index)~=-1)
    error('RocketAeroHeat:EulerEdge:Wall','需要至少三个连续的正半径机体壁面。');
end
if any(ds(index)<=0) || any(~isfinite(ds(index)))
    error('RocketAeroHeat:EulerEdge:Geometry','壁面线段长度无效。');
end
lengths=ds(index);s=cumsum(lengths)-.5*lengths;
tx=-dx(index)./lengths;tr=-dr(index)./lengths;
nx=tr;nr=-tx; % 外法向指向流体；在圆柱上为 +R。
rho=result.density(index,1);p=result.pressure(index,1);
ux=result.velocityX(index,1);ur=result.velocityR(index,1);
sampleX=g.centerX(index,1);sampleR=g.centerR(index,1);
x=g.wallX(index);radius=g.wallR(index);fromTrace=false;
if isfield(result,'spatialOrder') && result.spatialOrder==2
    if ~isfield(result.wall,'trace')
        error('RocketAeroHeat:EulerEdge:Trace','线性重构解缺少壁面迹，不能静默改用首层单元。');
    end
    trace=result.wall.trace;
    for key={'density','pressure','velocityX','velocityR','x','radius','fractionFromTail'}
        validateattributes(trace.(key{1}),{'numeric'},{'real','finite','size',[g.nx,1]});
    end
    if any(trace.density<=0) || any(trace.pressure<=0) || any(trace.fractionFromTail<0 | trace.fractionFromTail>1)
        error('RocketAeroHeat:EulerEdge:Trace','壁面重构迹的密度、压力或面内位置无效。');
    end
    rho=trace.density(index);p=trace.pressure(index);ux=trace.velocityX(index);ur=trace.velocityR(index);
    x=trace.x(index);radius=trace.radius(index);sampleX=x;sampleR=radius;fromTrace=true;
    % r 加权面位置不是几何中点。鼻端首段应在 2/3 弧长附近，不能
    % 把重构速度放在 1/2 弧长位置再求梯度，否则会产生系统性偏差。
    s=cumsum(lengths)-trace.fractionFromTail(index).*lengths;
end
ut=ux.*tx+ur.*tr;un=ux.*nx+ur.*nr;
cp=fs.gamma*fs.R/(fs.gamma-1);T=p./(rho*fs.R);
speed=abs(ut);sound=sqrt(fs.gamma*fs.R*T);M=speed./sound;
H=cp*T+.5*speed.^2;sampleH=cp*T+.5*(ux.^2+ur.^2);
Hinf=cp*fs.temperature+.5*fs.speed^2;
p0=p.*(1+.5*(fs.gamma-1)*M.^2).^(fs.gamma/(fs.gamma-1));

%% 2. 在几何折角处断开导数；不跨过肩部差分，也不人为补驻点 U=0
turn=acosd(max(-1,min(1,tx(1:end-1).*tx(2:end)+tr(1:end-1).*tr(2:end))));
feature=false(n,1);j=find(turn>cornerAngle);feature(j)=true;feature(j+1)=true;
reverse=ut<0;
smooth=~feature & ~reverse;
dUds=localDerivative(s,ut,smooth);dpds=localDerivative(s,p,smooth);
wallPressure=result.wall.pressure(index);
% 未收敛阶段的反射壁面通量可能给出非正的压力估计；保留原值并标记，
% 不应阻止查看原始单元仍为正的场，也不能裁剪它来伪装可用状态。
validateattributes(wallPressure,{'numeric'},{'real','finite','size',[n,1]});

%% 3. 保存候选状态和独立诊断；不得用小残差替代速度/梯度精度验证
edge=struct('s',s,'x',x,'radius',radius, ...
    'wallIndex',index,'segmentLength',lengths,'tangent',[tx,tr],'normal',[nx,nr], ...
    'pressure',p,'temperature',T,'density',rho,'speed',speed,'signedSpeed',ut, ...
    'velocityX',ut.*tx,'velocityR',ut.*tr,'Mach',M,'speedOfSound',sound, ...
    'totalEnthalpy',H,'totalTemperature',H/cp,'totalPressure',p0, ...
    'dynamicViscosity',1.458e-6*T.^1.5./(T+110.4), ...
    'velocityGradient',dUds,'pressureGradient',dpds);
edge.sample=struct('velocityX',ux,'velocityR',ur,'normalVelocity',un, ...
    'totalEnthalpy',sampleH,'centerX',sampleX,'centerR',sampleR);
edge.sample.fromReconstructedTrace=fromTrace;
edge.wallFluxPressure=wallPressure;
edge.diagnostics.wallPressureRelativeDifference=(wallPressure-p)./p;
edge.diagnostics.normalVelocityFraction=abs(un)./max(hypot(ux,ur),eps(fs.speed));
edge.diagnostics.projectionEnthalpyLoss=.5*un.^2/Hinf;
edge.diagnostics.sampleTotalEnthalpyRelativeError=sampleH/Hinf-1;
edge.diagnostics.projectedTotalEnthalpyRelativeError=H/Hinf-1;
edge.diagnostics.entropyChange=cp*log(T/fs.temperature)-fs.R*log(p/fs.pressure);
edge.diagnostics.totalPressureRatio=p0/fs.totalPressure;
edge.flags.feature=feature;edge.flags.reverseFlow=reverse;
edge.flags.wallFluxPressureNonpositive=wallPressure<=0;
edge.flags.gradientAvailable=isfinite(dUds) & isfinite(dpds);
edge.flags.cornerAngleDeg=cornerAngle;
edge.freestream=fs;edge.contour=g.contour;
edge.sourceMethod=result.method;
edge.thermodynamicModel='caloricallyPerfect';
edge.method='第一层 Euler 单元原始状态的一阶壁面迹近似；几何切向投影；定比热局部总状态';
if fromTrace
    edge.method='受限线性重构的壁面内侧迹；r 加权面位置；几何切向投影；定比热局部总状态';
end
edge.gradientMethod='沿鼻端向尾部弧长的非均匀三点二次差分；折角/逆流处分段，端点单边差分';
edge.scope='仅零攻角回转体无粘近壁候选；不是黏性边界层厚度处的实测状态';
edge.inputConverged=result.converged;edge.inputBoundaryClear=result.boundary.clear;
edge.algebraicallyValid=isfinite(H) & isfinite(p0) & p0>0;
edge.available=false;edge.gridAccuracyVerified=false;
edge.status='外缘速度、熵及梯度尚未通过网格精度验证；只作诊断，未接入热流';
end

function derivative=localDerivative(s,y,valid)
% 模块作用：在每个连续有效段内计算导数，避免跨 NaN/折角制造假梯度。
derivative=nan(size(y));starts=find(valid & [true;~valid(1:end-1)]);
ends=find(valid & [~valid(2:end);true]);
for k=1:numel(starts)
    first=starts(k);last=ends(k);
    if last-first<2,continue;end
    for i=first:last
        center=min(max(i,first+1),last-1);ids=center-1:center+1;
        delta=s(ids)-s(i);scale=max(abs(delta));z=delta/scale;
        % 先去掉常量以减少大静压上的差分相消误差。
        fit=[ones(3,1),z,z.^2]\(y(ids)-y(i));
        derivative(i)=fit(2)/scale;
    end
end
end
