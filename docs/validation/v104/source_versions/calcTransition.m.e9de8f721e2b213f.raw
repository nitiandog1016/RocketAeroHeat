function transition = calcTransition(edgeReynolds,eligibleFace,cfg,path,resultMask)
%CALCTRANSITION 沿表面准流线的 Reynolds 数经验转捩判据。
% 输入 Re_s=rho_e*U_e*s/mu_e，以及流动状态可传递的面元掩码。
% Re_s>=Re_crit 首次判湍流；若提供三维路径，下游继承已发生的转捩。
% 默认 Re_crit=5e5。这是固定阈值工程筛选，不能预测真实转捩位置、
% 湍流分离或高超声速第二模态；不提供路径时保留逐面元旧接口。
% 可选 resultMask 仅筛选最终可输出热流的面元，不切断其余有效外流
% 面元上传递的已转捩状态。例如热物性空白处可传递湍流标签，但不输出热流。

%% 1. 参数和数据校验
if nargin<3, cfg=struct; end
if nargin<4, path=[]; end
if nargin<5, resultMask=eligibleFace; end
critical=5e5;
if isfield(cfg,'heating') && isfield(cfg.heating,'transitionRe')
    critical=cfg.heating.transitionRe;
end
validateattributes(critical,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(edgeReynolds,{'numeric'},{'real','column'});
validateattributes(eligibleFace,{'logical','numeric'}, ...
    {'binary','size',size(edgeReynolds)});
eligibleFace=logical(eligibleFace);
validateattributes(resultMask,{'logical','numeric'}, ...
    {'binary','size',size(edgeReynolds)});
resultMask=logical(resultMask);
if any(resultMask & ~eligibleFace)
    error('RocketAeroHeat:Transition:Mask', ...
        '热流输出面元必须属于可传递流动状态的面元。');
end
if any(~isfinite(edgeReynolds(eligibleFace)) | edgeReynolds(eligibleFace)<=0)
    error('RocketAeroHeat:Transition:Reynolds', ...
        '可计算面元的局部 Reynolds 数必须为正有限值。');
end

%% 2. 先判局部阈值，再按有向路径传递已有湍流状态
stateTurbulent=eligibleFace & edgeReynolds>=critical;
stateLocalThreshold=stateTurbulent;
pathUsed=~isempty(path);
throughHidden=false(size(edgeReynolds));
if pathUsed
    required={'predecessorFace','pathCost','connectedFace'};
    for j=1:numel(required)
        if ~isfield(path,required{j}) || numel(path.(required{j}))~=numel(edgeReynolds)
            error('RocketAeroHeat:Transition:Path', ...
                '流线路径缺少 %s 或面元数不匹配。',required{j});
        end
    end
    predecessor=double(path.predecessorFace(:));
    cost=double(path.pathCost(:));
    connected=logical(path.connectedFace(:));
    if any(predecessor<0 | predecessor>numel(edgeReynolds) | ...
            predecessor~=floor(predecessor))
        error('RocketAeroHeat:Transition:Path','前驱面元编号无效。');
    end
    [~,order]=sort(cost(eligibleFace & connected));
    ids=find(eligibleFace & connected); ids=ids(order);
    for k=1:numel(ids)
        id=ids(k); parent=predecessor(id);
        if parent>0 && eligibleFace(parent) && ...
                isfinite(cost(parent)) && cost(parent)<cost(id) && ...
                stateTurbulent(parent)
            stateTurbulent(id)=true;
            throughHidden(id)=throughHidden(parent) || ~resultMask(parent);
        end
    end
end
stateLaminar=eligibleFace & ~stateTurbulent;
stateOnset=stateTurbulent;
if pathUsed
    children=find(stateTurbulent & connected & predecessor>0);
    stateOnset(children)=~stateTurbulent(predecessor(children));
end
% 输出掩码只影响面元结果，不影响上述有向路径的状态传播。
turbulent=resultMask & stateTurbulent;
laminar=resultMask & stateLaminar;
onset=resultMask & stateOnset;
localThreshold=resultMask & stateLocalThreshold;
intermittency=nan(size(edgeReynolds));
intermittency(laminar)=0; intermittency(turbulent)=1;
transition=struct;
if pathUsed
    transition.method='外缘 Reynolds 数固定阈值 + 沿三维准流线保持已发生的转捩';
else
    transition.method='局部外缘 Reynolds 数固定阈值经验判据';
end
transition.scope='附着流初步筛选；未使用扰动谱、粗糙度、可压缩稳定性或 e^N 分析';
transition.continuityMethod='先在有效外流路径传递转捩状态，再筛选热流输出';
transition.continuityScope=['热流空白只允许传递已有的层/湍流标签；' ...
    '不补热物性或动量厚度，不跨越外流、连续介质或路径失效；' ...
    '不表示已求出分离、再附或重新层流化'];
transition.criticalReynolds=critical;
transition.edgeReynolds=edgeReynolds;
transition.laminarFace=laminar;
transition.turbulentFace=turbulent;
transition.intermittency=intermittency;
transition.onsetFace=onset;
transition.localThresholdFace=localThreshold;
transition.stateEligibleFace=eligibleFace;
transition.resultEligibleFace=resultMask;
transition.stateTurbulentFace=stateTurbulent;
transition.stateOnsetFace=stateOnset;
transition.turbulentCarrierFace=stateTurbulent & ~resultMask;
transition.inheritedThroughUnavailableFace=turbulent & throughHidden;
transition.pathUsed=pathUsed;
transition.unresolvedFace=~resultMask;
transition.summary=struct('laminarFaces',nnz(laminar), ...
    'turbulentFaces',nnz(turbulent),'unresolvedFaces',nnz(~resultMask), ...
    'onsetFaces',nnz(onset), ...
    'inheritedTurbulentFaces',nnz(turbulent & ~localThreshold), ...
    'turbulentCarrierFaces',nnz(transition.turbulentCarrierFace), ...
    'inheritedThroughUnavailableFaces', ...
    nnz(transition.inheritedThroughUnavailableFace), ...
    'onsetsOutsideHeatDomain',nnz(stateOnset & ~resultMask));
end
