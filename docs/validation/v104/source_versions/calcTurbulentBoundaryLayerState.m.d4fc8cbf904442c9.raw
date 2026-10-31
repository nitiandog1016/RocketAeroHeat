function result = calcTurbulentBoundaryLayerState(mesh,nuReference,transition,cfg)
%CALCTURBULENTBOUNDARYLAYERSTATE 沿三维准流线推进准二维湍流积分模型。
% 模块作用：默认 Head 卷吸 + von Karman 动量积分，H 随流动变化；
% Cf 用 Ludwieg-Tillmann，参考运动黏度由 Eckert 参考温度提供。
% 层流交接/缺解后的局部平板启动采用 H=1.4，属于工程初值。
% H>=2.4 视为疑似准二维湍流分离；没有再附模型，整条下游无解。
% Re_theta 域外仍保留显式平板回退；数值失败单独标记并阻断下游。
% 在相关式内侧 300~600 与 8e5~1e6 的带内，用 smoothstep 权重将
% theta/H/Cf 连续衔接到平板状态。带宽是工程数值选择，不是转捩判据；
% 域外不外推 LT。平板回退保留平板 theta/H，避免下游状态硬切换。
% cfg.heating.turbulentModel='fixedShape' 可复现旧版固定 H 模型。

%% 1. 选择模型；复用旧模块的输入校验和可审计平板基准
if nargin<4, cfg=struct; end
model='head';
if isfield(cfg,'heating') && isfield(cfg.heating,'turbulentModel')
    model=char(cfg.heating.turbulentModel);
end
if ~ismember(model,{'head','fixedShape'})
    error('RocketAeroHeat:TurbulentBL:Model','湍流模型须为 head 或 fixedShape。');
end
result=calcFixedShapeTurbulentBoundaryLayerState(mesh,nuReference,transition);
nf=size(mesh.F,1);
result.model=model;
result.possibleSeparationFace=false(nf,1);
result.separationUnresolvedFace=false(nf,1);
result.numericalFailureFace=false(nf,1);
result.numericalUnresolvedFace=false(nf,1);
result.separationThreshold=2.4;
result.shapeFactor=1.4*ones(nf,1);
result.shapeFactor(~isfinite(result.momentumThickness))=NaN;
result.summary.separationOnsetFaces=0;
result.summary.separationUnresolvedFaces=0;
result.summary.numericalFailureFaces=0;
result.summary.numericalUnresolvedFaces=0;
result.summary.blendedFaces=0;
if strcmp(model,'fixedShape'), return; end

%% 2. 沿前驱图按上游到下游顺序推进，分离状态穿过热物性空白传递
p=mesh.surfaceStreamlines; U=double(mesh.boundaryLayerEdge.speed(:));
nu=double(nuReference(:)); parent=double(p.predecessorFace(:));
cost=double(p.pathCost(:)); connected=logical(p.connectedFace(:));
turbulent=logical(transition.turbulentFace(:));
laminar=logical(transition.laminarFace(:));
theta=nan(nf,1); H=theta; Cf=theta; gradient=theta; beta=theta;
source=zeros(nf,1,'uint8');
weight=nan(nf,1);
separated=false(nf,1); failure=separated;
onset=separated; failOnset=separated;
flatCf=result.flatPlateSkinFrictionCoefficient;
ids=find(connected); [~,order]=sort(cost(ids)); ids=ids(order);
for k=1:numel(ids)
    id=ids(k); up=parent(id);
    if up>0
        if ~connected(up) || cost(up)>=cost(id)
            error('RocketAeroHeat:TurbulentBL:Path','前驱必须是成本更小的有效上游面元。');
        end
        separated(id)=separated(up); failure(id)=failure(up);
    end
    % 即使此面元不输出热流，也必须保持上游分离/数值失效状态。
    if separated(id), source(id)=5; continue; end
    if failure(id), source(id)=6; continue; end
    if ~turbulent(id) || ~isfinite(flatCf(id)) || ...
            ~isfinite(U(id)) || U(id)<=0 || ~isfinite(nu(id)) || nu(id)<=0
        continue;
    end
    th0=NaN; H0=1.4; kind=uint8(0);
    Re=U(id)*p.distance(id)/nu(id);
    flatTheta=0.037*p.distance(id)/Re^0.2;
    if up>0 && isfinite(U(up)) && U(up)>0 && isfinite(nu(up)) && nu(up)>0
        if turbulent(up) && isfinite(theta(up)) && isfinite(H(up))
            th0=theta(up); H0=H(up); kind=uint8(2);
        elseif laminar(up) && isfinite(mesh.boundaryLayer.momentumThickness(up))
            th0=mesh.boundaryLayer.momentumThickness(up); kind=uint8(1);
        end
    end
    domainFallback=false;
    if isfinite(th0) && th0>0 && isfinite(p.stepLength(id)) && p.stepLength(id)>0
        segment=marchHeadBoundaryLayerSegment(th0,H0,U(up),U(id), ...
            nu(up),nu(id),p.stepLength(id));
        if segment.status==0
            weight(id)=domainWeight(segment.reynoldsRange(1),segment.reynoldsRange(2));
            w=weight(id);
            theta(id)=w*segment.theta+(1-w)*flatTheta;
            H(id)=w*segment.shapeFactor+(1-w)*1.4;
            Cf(id)=w*segment.Cf+(1-w)*flatCf(id); source(id)=kind;
            if w<1, source(id)=7; end
            gradient(id)=(U(id)-U(up))/p.stepLength(id);
            beta(id)=-2*H(id)*theta(id)*gradient(id)/(Cf(id)*U(id));
            continue;
        elseif segment.status==1
            onset(id)=true; separated(id)=true; source(id)=5; continue;
        elseif segment.status==3
            failOnset(id)=true; failure(id)=true; source(id)=6; continue;
        end
        domainFallback=true;
    end
    % 没有连续上游解时才局部启动；相关式域外显式回退，不冒充积分。
    candidate=flatTheta;
    reTheta=U(id)*candidate/nu(id);
    if ~domainFallback && isfinite(reTheta) && reTheta>=300 && reTheta<=1e6
        theta(id)=candidate; H(id)=1.4;
        weight(id)=domainWeight(reTheta,reTheta);
        Cf(id)=weight(id)*0.246*10^(-0.678*H(id))*reTheta^(-0.268)+ ...
            (1-weight(id))*flatCf(id); source(id)=3;
    else
        Cf(id)=flatCf(id); source(id)=4;
        theta(id)=candidate; H(id)=1.4; weight(id)=0;
    end
end

%% 3. 保存逐面元状态、算法和范围；分离与积分失败不混为一类
result.method='Head 卷吸（H1 连续化）+ von Karman/Ludwieg-Tillmann（可变 H，域边平板衔接）';
result.scope=['沿三维准流线的准二维工程闭合，使用 Eckert 参考黏度；' ...
    '300<=Re_theta<=1e6，H>=2.4 仅为疑似分离；无再附模型，分离及下游无解；' ...
    'Re_theta 内侧 300~600、8e5~1e6 与平板连续衔接（工程带宽），' ...
    '衔接状态不冒充纯 Head 解；局部启动/域外平板 H=1.4；非完整可压缩边界层，' ...
    '不预测三维横流、真实分离线或激波干扰'];
result.shapeFactor=H;
result.momentumThickness=theta;
result.momentumThicknessReynolds=U.*theta./nu;
result.skinFrictionCoefficient=Cf;
result.edgeVelocityGradient=gradient;
result.adversePressureGradientParameter=beta;
result.sourceID=source;
result.integralWeight=weight;
result.blendingReynoldsBands=[300 600;8e5 1e6];
result.sourceNames={'非湍流或无解','从层流转捩接续','Head 湍流逐段积分', ...
    '局部湍流平板启动','LT 域外平板回退','疑似湍流分离及下游','Head 积分失败及下游', ...
    'LT 域边与平板连续衔接'};
result.possibleSeparationFace=onset;
result.separationUnresolvedFace=separated;
result.numericalFailureFace=failOnset;
result.numericalUnresolvedFace=failure;
result.validFace=turbulent & isfinite(Cf) & Cf>0 & ~separated & ~failure;
result.fallbackFace=source==4;
result.summary=struct('validFaces',nnz(result.validFace), ...
    'laminarHandoffFaces',nnz(source==1),'integratedFaces',nnz(source==2), ...
    'localStartupFaces',nnz(source==3),'flatPlateFallbackFaces',nnz(source==4), ...
    'adverseGradientFaces',nnz(beta>0),'separationOnsetFaces',nnz(onset), ...
    'separationUnresolvedFaces',nnz(separated),'numericalFailureFaces',nnz(failOnset), ...
    'numericalUnresolvedFaces',nnz(failure),'blendedFaces',nnz(source==7));
end

function w=domainWeight(low,high)
% 两端取三次 smoothstep，权重及一阶导数在衔接端点连续。
a=max(0,min(1,(low-300)/300));
b=max(0,min(1,(1e6-high)/2e5));
w=a^2*(3-2*a)*b^2*(3-2*b);
end
