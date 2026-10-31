function mesh = assessContinuumRegime(mesh,freestream,cfg)
%ASSESSCONTINUUMREGIME 筛选连续介质气动热公式的尺度适用性。
% 模块作用：用自由来流和外缘的平均自由程除以同一特征长度，
% 给出逐面元 Knudsen 数；有效状态达到阈值即不再套用当前无滑移
% 激波/边界层/热流工程模型。未解热力状态单独阻断并沿准流线传播，
% 域外流向占位不得标成真实稀薄流，也不能用未知状态放行热流。
% 默认 Kn<0.01 是保守的连续介质筛选，不是局部梯度 Knudsen 数，
% 也不提供滑移边界、过渡流或自由分子热流解。

%% 1. 校验外缘状态和几何特征尺度
if nargin<3, cfg=struct; end
if ~isfield(mesh,'boundaryLayerEdge') || ~isfield(mesh,'F')
    error('RocketAeroHeat:Continuum:Input', ...
        '需要当前网格及边界层外缘压力、温度和黏度。');
end
nf=size(mesh.F,1);
e=mesh.boundaryLayerEdge;
needed={'pressure','temperature','dynamicViscosity','algebraicallyValid'};
for j=1:numel(needed)
    if ~isfield(e,needed{j}) || numel(e.(needed{j}))~=nf
        error('RocketAeroHeat:Continuum:Edge', ...
            '外缘状态缺少逐面元 %s。',needed{j});
    end
end
limit=0.01;
if isfield(cfg,'heating') && isfield(cfg.heating,'maxContinuumKnudsen')
    limit=cfg.heating.maxContinuumKnudsen;
end
validateattributes(limit,{'numeric'}, ...
    {'real','finite','scalar','>',0,'<=',0.01});
if isfield(cfg,'heating') && isfield(cfg.heating,'continuumLength')
    lengthScale=cfg.heating.continuumLength;
    lengthSource='用户指定的连续介质判据特征长度';
elseif isfield(mesh,'stagnationHeating') && ...
        isfinite(mesh.stagnationHeating.radius) && ...
        mesh.stagnationHeating.radius>0
    lengthScale=mesh.stagnationHeating.radius;
    lengthSource='已确认或经多尺度检查的物理鼻端半径';
elseif isfield(mesh,'stagnation') && ...
        isfield(mesh.stagnation,'hasPhysicalRadius') && ...
        mesh.stagnation.hasPhysicalRadius
    lengthScale=mesh.stagnation.physicalRadius;
    lengthSource='驻点模块给定的物理鼻端半径';
else
    if ~isfield(mesh,'V') || isempty(mesh.V)
        error('RocketAeroHeat:Continuum:Geometry', ...
            '没有物理鼻端半径或可用于估计横向尺度的网格顶点。');
    end
    span=max(mesh.V,[],1)-min(mesh.V,[],1);
    lengthScale=0.5*max(span(2:3));
    lengthSource='最大横向半径的替代尺度；不是已验证鼻端半径';
end
validateattributes(lengthScale,{'numeric'}, ...
    {'real','finite','scalar','positive'});
for name={'pressure','temperature','R'}
    if ~isfield(freestream,name{1})
        error('RocketAeroHeat:Continuum:Freestream', ...
            '自由来流缺少 %s。',name{1});
    end
end

%% 2. 计算自由来流和外缘 Kn；较大的值作为逐面元筛选指标
muInf=NaN;
if isfield(freestream,'dynamicViscosity')
    muInf=freestream.dynamicViscosity;
end
if ~isfinite(muInf) && isfinite(freestream.temperature) && ...
        freestream.temperature>0
    T=freestream.temperature;
    muInf=1.458e-6*T^(3/2)/(T+110.4);
end
lambdaInf=airMeanFreePath(freestream.pressure, ...
    freestream.temperature,muInf,freestream.R);
lambdaEdge=airMeanFreePath(e.pressure(:),e.temperature(:), ...
    e.dynamicViscosity(:),freestream.R);
knInf=lambdaInf/lengthScale;
knEdge=lambdaEdge/lengthScale;
valid=logical(e.algebraicallyValid(:));
stateValid=valid;
if isfield(e,'thermalStateValidFace')
    if numel(e.thermalStateValidFace)~=nf
        error('RocketAeroHeat:Continuum:StateSize','外缘物理场有效标志与面元不一致。');
    end
    stateValid=stateValid & logical(e.thermalStateValidFace(:));
end
% 兼容旧保存快照：旧占位也不能因没有新字段被认作真实低温外缘。
if isfield(e,'thermalExpansionOutOfRangeFace')
    if numel(e.thermalExpansionOutOfRangeFace)~=nf
        error('RocketAeroHeat:Continuum:StateSize','域外膨胀标志与面元不一致。');
    end
    stateValid=stateValid & ~logical(e.thermalExpansionOutOfRangeFace(:));
end
diagnosticKn=max(knEdge,knInf);
diagnosticKn(~valid)=NaN;
lambdaEdge(~stateValid)=NaN;knEdge(~stateValid)=NaN;
kn=max(knEdge,knInf);
kn(~valid)=NaN;
kn(~stateValid)=NaN;
directStateUnresolved=valid & (~stateValid | ~isfinite(kn));
directKnUnresolved=valid & stateValid & isfinite(kn) & kn>=limit;
stateBlocked=directStateUnresolved;knBlocked=directKnUnresolved;

%% 3. 上游失效后，不把同一条准流线下游重新当作附着连续介质解
directUnresolved=directStateUnresolved | directKnUnresolved;
if isfield(mesh,'surfaceStreamlines')
    path=mesh.surfaceStreamlines;
    if ~isfield(path,'predecessorFace') || ...
            ~isfield(path,'pathCost') || ...
            ~isfield(path,'connectedFace') || ...
            ~isfield(path,'velocitySource') || ...
            numel(path.predecessorFace)~=nf || ...
            numel(path.pathCost)~=nf || numel(path.connectedFace)~=nf || ...
            ~strcmp(path.velocitySource,e.method)
        error('RocketAeroHeat:Continuum:Path', ...
            '准流线与当前外缘状态不一致。');
    end
    connected=logical(path.connectedFace(:));
    parent=double(path.predecessorFace(:));
    cost=double(path.pathCost(:));
    if any(parent<0 | parent>nf | parent~=floor(parent)) || ...
            any(~isfinite(cost(connected)))
        error('RocketAeroHeat:Continuum:Path', ...
            '准流线前驱或路径成本无效。');
    end
    ids=find(connected);
    [~,order]=sort(cost(ids)); ids=ids(order);
    for j=1:numel(ids)
        id=ids(j); upstream=parent(id);
        if upstream==0, continue; end
        if ~connected(upstream) || cost(upstream)>=cost(id)
            error('RocketAeroHeat:Continuum:Path', ...
                '准流线前驱必须是更上游的有效面元。');
        end
        % 两类风险各自传播，不能把域外占位产生的巨大Kn当作真实稀薄区。
        stateBlocked(id)=stateBlocked(id) || stateBlocked(upstream);
        knBlocked(id)=knBlocked(id) || knBlocked(upstream);
    end
end
unresolved=stateBlocked | knBlocked;

%% 4. 保存逐面元区域、方法和适用范围供热流与 GUI 使用
regime=zeros(nf,1,'uint8');
regime(stateValid & kn<0.01)=1;
regime(stateValid & kn>=0.01 & kn<0.1)=2;
regime(stateValid & kn>=0.1 & kn<10)=3;
regime(stateValid & kn>=10)=4;
r=struct;
r.method='黏度动理论平均自由程 / 物理鼻端或指定特征长度 Kn 筛选';
r.scope=['只有有效外缘状态才有Kn=max(Kn_{自由来流},Kn_{外缘})；默认 Kn<0.01 才使用当前' ...
    '无滑移连续介质热流模型；这不是梯度 Kn、滑移流或 DSMC 解；' ...
    '域外热力占位及下游单独标为外流未解，仍阻断热流；' ...
    '若用横向半径替代鼻端半径，阈值仅供初步筛选'];
r.characteristicLength=lengthScale;
r.lengthSource=lengthSource;
r.maxContinuumKnudsen=limit;
r.freestreamMeanFreePath=lambdaInf;
r.freestreamKnudsen=knInf;
r.edgeMeanFreePath=lambdaEdge;
r.edgeKnudsen=knEdge;
r.knudsen=kn;
r.diagnosticKnudsen=diagnosticKn;
r.diagnosticKnudsenScope='包含仅供流向诊断的域外定比热占位；不得作为真实流态或热流有效性证据';
r.thermalStateValidFace=stateValid;
r.regimeID=regime;
r.regimeNames={'无有效外缘','连续介质','滑移流','过渡流','自由分子流'};
r.directUnresolvedFace=directUnresolved;
r.downstreamUnresolvedFace=unresolved & ~directUnresolved;
r.directKnudsenUnresolvedFace=directKnUnresolved;
r.downstreamKnudsenUnresolvedFace=knBlocked & ~directKnUnresolved;
r.knudsenUnresolvedFace=knBlocked;
r.directStateUnresolvedFace=directStateUnresolved;
r.downstreamStateUnresolvedFace=stateBlocked & ~directStateUnresolved;
r.stateUnresolvedFace=stateBlocked;
r.unresolvedFace=unresolved;
% 下游可重新呈现较小的局部 Kn，但上游已离开连续介质适用域。
% 单独给出显示编号，避免绿色“连续介质”误导为可计算热流。
r.applicabilityID=regime;
r.applicabilityID(r.downstreamKnudsenUnresolvedFace)=5;
r.applicabilityID(r.downstreamStateUnresolvedFace)=7;
r.applicabilityID(r.directStateUnresolvedFace)=6;
r.applicabilityNames=[r.regimeNames,{'上游有效外缘的Kn限制', ...
    '外缘热力状态未解/占位','上游外缘热力状态未解/占位'}];
r.summary=struct('continuumFaces',nnz(valid & ~unresolved), ...
    'directUnresolvedFaces',nnz(directUnresolved), ...
    'downstreamUnresolvedFaces',nnz(r.downstreamUnresolvedFace), ...
    'unknownKnudsenFaces',nnz(valid & ~isfinite(kn)), ...
    'directKnudsenUnresolvedFaces',nnz(directKnUnresolved), ...
    'downstreamKnudsenUnresolvedFaces',nnz(r.downstreamKnudsenUnresolvedFace), ...
    'directStateUnresolvedFaces',nnz(directStateUnresolved), ...
    'downstreamStateUnresolvedFaces',nnz(r.downstreamStateUnresolvedFace));
mesh.continuumRegime=r;
end
