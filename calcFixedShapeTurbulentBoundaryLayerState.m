function result = calcFixedShapeTurbulentBoundaryLayerState(mesh,nuReference,transition)
%CALCTURBULENTBOUNDARYLAYERSTATE 沿三维准流线积分附着湍流动量厚度。
% 模块作用：从层流转捩面元的 theta 接续，按 von Karman 动量积分式
% d(theta*Ue^(H+2))/ds=(Cf/2)*Ue^(H+2) 推进；局部摩擦系数用
% Ludwieg-Tillmann 关联式 Cf=0.246*10^(-0.678H)*Re_theta^(-0.268)。
% 固定 H=1.4，参考运动黏度与热流模块的 Eckert 湍流参考温度一致。
% 缺少有效上游边界层时以局部湍流平板 theta 启动；相关式超出
% 300<=Re_theta<=1e6 时明确退回原局部平板 Cf 并标记来源。
% 这是沿离散准流线的准二维、附着湍流工程近似，不求解三维横流、
% 湍流分离、再附、强激波干扰或真实可压缩边界层全剖面。

%% 1. 检查外缘、路径、层流交接值与转捩面元
if ~isfield(mesh,'boundaryLayerEdge') || ...
        ~isfield(mesh,'surfaceStreamlines') || ...
        ~isfield(mesh,'boundaryLayer')
    error('RocketAeroHeat:TurbulentBL:Input', ...
        '需要外缘状态、三维准流线和层流动量厚度。');
end
nf=size(mesh.F,1);
e=mesh.boundaryLayerEdge;
path=mesh.surfaceStreamlines;
laminar=mesh.boundaryLayer;
nuReference=double(nuReference(:));
required={'predecessorFace','stepLength','pathCost','distance', ...
    'connectedFace','velocitySource'};
for j=1:numel(required)
    if ~isfield(path,required{j})
        error('RocketAeroHeat:TurbulentBL:Path', ...
            '准流线缺少 %s。',required{j});
    end
end
if numel(nuReference)~=nf || numel(e.speed)~=nf || ...
        numel(path.predecessorFace)~=nf || ...
        numel(path.stepLength)~=nf || numel(path.pathCost)~=nf || ...
        numel(path.distance)~=nf || numel(path.connectedFace)~=nf || ...
        numel(laminar.momentumThickness)~=nf || ...
        numel(transition.turbulentFace)~=nf || ...
        numel(transition.laminarFace)~=nf || ...
        ~strcmp(path.velocitySource,e.method)
    error('RocketAeroHeat:TurbulentBL:Size', ...
        '逐面元外缘、层流、转捩或准流线与当前网格不一致。');
end
parent=double(path.predecessorFace(:));
cost=double(path.pathCost(:));
distance=double(path.distance(:));
stepLength=double(path.stepLength(:));
connected=logical(path.connectedFace(:));
turbulent=logical(transition.turbulentFace(:));
laminarFace=logical(transition.laminarFace(:));
U=double(e.speed(:));
if any(parent<0 | parent>nf | parent~=floor(parent)) || ...
        any(~isfinite(cost(connected)))
    error('RocketAeroHeat:TurbulentBL:Path', ...
        '准流线前驱编号或路径成本无效。');
end

%% 2. 给出局部平板对照值；只在湍流面元沿前驱树积分
H=1.4;
power=H+2;
coefficient=0.246*10^(-0.678*H);
minimumReTheta=300;
maximumReTheta=1e6;
ReDistance=U.*distance./nuReference;
flatCf=nan(nf,1);
flatValid=turbulent & connected & isfinite(ReDistance) & ...
    ReDistance>0;
flatCf(flatValid)=0.0592./ReDistance(flatValid).^0.2;
theta=nan(nf,1);
ReTheta=theta;
Cf=theta;
gradient=theta;
beta=theta;
source=zeros(nf,1,'uint8');
ids=find(turbulent & connected);
[~,order]=sort(cost(ids));
ids=ids(order);
for j=1:numel(ids)
    id=ids(j);
    if ~isfinite(U(id)) || U(id)<=0 || ...
            ~isfinite(nuReference(id)) || nuReference(id)<=0 || ...
            ~isfinite(flatCf(id))
        continue;
    end
    upstream=parent(id);
    thetaParent=NaN;
    cfParent=NaN;
    if upstream>0 && connected(upstream) && ...
            cost(upstream)<cost(id) && ...
            isfinite(U(upstream)) && U(upstream)>0 && ...
            isfinite(nuReference(upstream)) && ...
            nuReference(upstream)>0
        if turbulent(upstream) && isfinite(theta(upstream)) && ...
                isfinite(Cf(upstream))
            thetaParent=theta(upstream);
            cfParent=Cf(upstream);
            source(id)=2;
        elseif laminarFace(upstream) && ...
                isfinite(laminar.momentumThickness(upstream)) && ...
                laminar.momentumThickness(upstream)>0
            thetaParent=laminar.momentumThickness(upstream);
            reParent=U(upstream)*thetaParent/nuReference(upstream);
            if reParent>=minimumReTheta && reParent<=maximumReTheta
                cfParent=coefficient*reParent^(-0.268);
                source(id)=1;
            end
        end
    end
    if isfinite(thetaParent) && isfinite(cfParent) && ...
            isfinite(stepLength(id)) && stepLength(id)>0
        ds=stepLength(id);
        ratio=U(upstream)/U(id);
        transported=thetaParent*ratio^power;
        candidate=transported+ds*cfParent*(ratio^power+1)/4;
        for iteration=1:8
            re=U(id)*candidate/nuReference(id);
            if ~(isfinite(re) && re>=minimumReTheta && ...
                    re<=maximumReTheta), break; end
            localCf=coefficient*re^(-0.268);
            next=transported+ds*(cfParent*ratio^power+localCf)/4;
            if abs(next-candidate)<=1e-9*max(candidate,eps)
                candidate=next;
                break;
            end
            candidate=next;
        end
        re=U(id)*candidate/nuReference(id);
        if isfinite(re) && re>=minimumReTheta && ...
                re<=maximumReTheta && isfinite(candidate) && candidate>0
            theta(id)=candidate;
            ReTheta(id)=re;
            Cf(id)=coefficient*re^(-0.268);
            gradient(id)=(U(id)-U(upstream))/ds;
            beta(id)=-2*H*theta(id)*gradient(id)/(Cf(id)*U(id));
            continue;
        end
    end
    % 上游状态不连续时，以湍流平板 theta 为本地初值；明确记录来源。
    if isfinite(distance(id)) && distance(id)>0
        trialTheta=0.037*distance(id)/ReDistance(id)^0.2;
        trialRe=U(id)*trialTheta/nuReference(id);
        if isfinite(trialRe) && trialRe>=minimumReTheta && ...
                trialRe<=maximumReTheta
            theta(id)=trialTheta;
            ReTheta(id)=trialRe;
            Cf(id)=coefficient*trialRe^(-0.268);
            source(id)=3;
            continue;
        end
    end
    % 相关式域外不外推 LT；保留可审计的旧平板值作为低阶回退。
    Cf(id)=flatCf(id);
    source(id)=4;
end

%% 3. 保存算法、来源与失效标志，供热流和 GUI 逐面元辨认
result=struct;
result.method=['沿三维准流线的 von Karman 湍流动量积分 + ' ...
    'Ludwieg-Tillmann Cf（固定 H=1.4）'];
result.scope=['准二维附着湍流工程模型；使用 Eckert 参考黏度，' ...
    '300<=Re_theta<=1e6；缺少上游解时局部启动，相关式域外退回' ...
    '局部平板 Cf；不预测三维横流、湍流分离或再附'];
result.shapeFactor=H;
result.correlationReynoldsRange=[minimumReTheta maximumReTheta];
result.momentumThickness=theta;
result.momentumThicknessReynolds=ReTheta;
result.skinFrictionCoefficient=Cf;
result.flatPlateSkinFrictionCoefficient=flatCf;
result.edgeVelocityGradient=gradient;
result.adversePressureGradientParameter=beta;
result.sourceID=source;
result.sourceNames={'非湍流或无解','从层流转捩接续','湍流逐段积分', ...
    '局部湍流平板启动','LT 域外平板回退'};
result.validFace=turbulent & isfinite(Cf) & Cf>0;
result.fallbackFace=source==4;
result.summary=struct('validFaces',nnz(result.validFace), ...
    'laminarHandoffFaces',nnz(source==1), ...
    'integratedFaces',nnz(source==2), ...
    'localStartupFaces',nnz(source==3), ...
    'flatPlateFallbackFaces',nnz(source==4), ...
    'adverseGradientFaces',nnz(beta>0));
end

