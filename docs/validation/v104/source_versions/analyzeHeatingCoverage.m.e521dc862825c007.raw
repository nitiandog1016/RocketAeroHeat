function mesh = analyzeHeatingCoverage(mesh)
%ANALYZEHEATINGCOVERAGE 解释面元热流有解范围及空白区的首要原因。
% 模块作用：只读取已经计算的几何、外流、准流线、边界层和热流标志，
% 对未解面元按固定优先级给出一个主要原因；不补值、不改变热流。
% 同一面元可能同时触发多个限制，flags 保留全部标志，reasonCode
% 只用于云图与面积统计。计算面积使用真实三角形面积，不用面元数量代替。

%% 1. 核对当前热流与网格，防止旧工况掩码混用
required={'heating','localGeometry','aerothermalRegion','boundaryLayerEdge'};
for j=1:numel(required)
    if ~isfield(mesh,required{j})
        error('RocketAeroHeat:HeatingCoverage:Input', ...
            '缺少 mesh.%s。',required{j});
    end
end
nf=size(mesh.F,1);
h=mesh.heating; g=mesh.localGeometry;
reg=mesh.aerothermalRegion; edge=mesh.boundaryLayerEdge;
area=double(g.area(:));
computed=logical(h.eligibleFace(:));
if numel(area)~=nf || numel(computed)~=nf || ...
        numel(edge.algebraicallyValid)~=nf || ...
        any(~isfinite(area) | area<=0) || ...
        ~isequal(logical(h.unresolvedFace(:)),~computed)
    error('RocketAeroHeat:HeatingCoverage:Size', ...
        '几何面积、外流或热流有效标志与面元不一致。');
end

%% 2. 收集可能重叠的物理/数值限制
flags=struct;
flags.endCap=logical(g.isAftCap(:)|g.isNoseCap(:)|g.isAftTipAdjacent(:));
flags.geometry=logical(g.topologyIssue(:)|g.featureFace(:));
flags.noseModel=logical(reg.isNoseTipAdjacent(:)|h.nearTipFace(:));
if isfield(mesh,'stagnation')
    flags.noseModel=flags.noseModel|logical(mesh.stagnation.candidateFace(:));
end
flags.externalFlow=~logical(edge.algebraicallyValid(:));
flags.rarefaction=false(nf,1);
if isfield(h,'continuumUnresolvedFace')
    flags.rarefaction=logical(h.continuumUnresolvedFace(:));
end
flags.upstreamExternalState=false(nf,1);
if isfield(mesh,'continuumRegime') && ...
        isfield(mesh.continuumRegime,'knudsenUnresolvedFace')
    c=mesh.continuumRegime;
    flags.rarefaction=logical(c.knudsenUnresolvedFace(:));
    flags.externalFlow=flags.externalFlow | logical(c.directStateUnresolvedFace(:));
    flags.upstreamExternalState=logical(c.downstreamStateUnresolvedFace(:));
end
flags.streamline=false(nf,1);
if isfield(mesh,'surfaceStreamlines')
    flags.streamline=~logical(mesh.surfaceStreamlines.connectedFace(:));
end
flags.laminarSeparation=false(nf,1);
flags.laminarCorrelation=false(nf,1);
flags.thermalProperty=false(nf,1);
flags.wallTemperature=false(nf,1);
% 显式选择的区域模型保留自身拒绝原因；旧准流线无效不能掩盖其物理域限制。
flags.regionalNoseRejected=false(nf,1);
regionalNoseCode=zeros(nf,1,'uint8');regionalNoseNames={};
if isfield(h,'noseModel') && strcmp(h.noseModel,'sphericalLocalSimilarity')
    nose=h.noseDistribution;
    flags.regionalNoseRejected=nose.edgeState.geometryFace(:) & nose.statusCode(:)>0;
    regionalNoseCode=nose.statusCode(:);regionalNoseNames=nose.statusNames(2:end);
end
if isfield(h,'wallTemperatureUnavailableFace')
    flags.wallTemperature=logical(h.wallTemperatureUnavailableFace(:));
end
flags.turbulentSeparation=false(nf,1);
flags.turbulentNumerical=false(nf,1);
if isfield(h,'turbulentSeparationUnresolvedFace')
    flags.turbulentSeparation=logical(h.turbulentSeparationUnresolvedFace(:));
    flags.turbulentNumerical=logical(h.turbulentNumericalUnresolvedFace(:));
end
if isfield(mesh,'boundaryLayer')
    b=mesh.boundaryLayer;
    flags.laminarSeparation=logical(b.possibleLaminarSeparationFace(:));
    flags.laminarCorrelation=logical(b.correlationOutOfRangeFace(:));
end
if isfield(h,'separationUnresolvedFace')
    flags.laminarSeparation=flags.laminarSeparation | ...
        logical(h.separationUnresolvedFace(:));
end
if isfield(h,'thermalPropertyUnavailableFace')
    flags.thermalProperty=logical(h.thermalPropertyUnavailableFace(:));
end
for key=fieldnames(flags)'
    if numel(flags.(key{1}))~=nf
        error('RocketAeroHeat:HeatingCoverage:Size', ...
            '%s 标志与网格面元数量不一致。',key{1});
    end
end

%% 3. 主要原因按可解释的优先级排他分类；最后保留“其他”而不猜测
% 代码 0=有解，1=端盖，2=特征/拓扑，3=鼻端模型，4=外流无解，
% 5=流线不连通，6=疑似层流分离，7=层流相关式范围外，
% 8=空气热物性范围外，9=其他，10=连续介质判据不适用，
% 11=疑似湍流分离及下游，12=Head 数值/形状域失效及下游，
% 13=壁温未知或超出热物性范围。
names={'已计算','端盖/尾端','特征折角或拓扑问题', ...
    '鼻端模型不适用','外流无有效解','准流线不连通', ...
    '疑似层流分离及下游','层流相关式范围外', ...
    '空气热物性模型范围外','其他未解限制', ...
    '连续介质假设不适用','疑似湍流分离及下游','Head 积分失败及下游', ...
    '壁温未知或超出热物性范围','上游外缘热力状态未解/域外占位影响'};
regionalOffset=numel(names)-1;
for j=1:numel(regionalNoseNames)
    names{end+1}=['所选球鼻模型：' regionalNoseNames{j}]; %#ok<AGROW>
end
code=zeros(nf,1,'uint8');
pending=~computed;
hit=pending & flags.regionalNoseRejected;
code(hit)=uint8(regionalOffset+double(regionalNoseCode(hit)));
pending(hit)=false;
masks={flags.endCap,flags.geometry,flags.noseModel, ...
    flags.rarefaction,flags.externalFlow,flags.upstreamExternalState,flags.streamline, ...
    flags.wallTemperature,flags.turbulentSeparation,flags.turbulentNumerical, ...
    flags.laminarSeparation,flags.laminarCorrelation,flags.thermalProperty};
priorityCode=[1 2 3 10 4 14 5 13 11 12 6 7 8];
for j=1:numel(masks)
    hit=pending & masks{j};
    code(hit)=uint8(priorityCode(j));
    pending(hit)=false;
end
code(pending)=9;
countByReason=accumarray(double(code)+1,1,[numel(names) 1])';
areaByReason=accumarray(double(code)+1,area,[numel(names) 1])';
totalArea=sum(area);

%% 4. 输出面积覆盖率、逐面元代码和算法适用范围
coverage=struct;
coverage.method='现有几何/外流/准流线/边界层有效标志的优先级诊断';
coverage.scope=['未解原因是工程诊断，不是新的热流模型；' ...
    '多个限制同时存在时只显示首要原因，全部原始标志保存在 flags'];
coverage.reasonCode=code;
coverage.reasonNames=names;
coverage.flags=flags;
coverage.countByReason=countByReason;
coverage.areaByReason=areaByReason;
coverage.summary=struct('totalFaces',nf, ...
    'computedFaces',nnz(computed), ...
    'unresolvedFaces',nnz(~computed), ...
    'totalArea',totalArea, ...
    'computedArea',areaByReason(1), ...
    'unresolvedArea',totalArea-areaByReason(1), ...
    'computedAreaFraction',areaByReason(1)/totalArea, ...
    'unresolvedAreaFraction',1-areaByReason(1)/totalArea);
mesh.heatingCoverage=coverage;
end
