function mesh = classifyAerothermalRegion(mesh,freestream,cfg)
%CLASSIFYAEROTHERMALREGION 根据最终网格与来流方向进行几何区域分类。
%
% 用途：为后续压力和气动热模块提供迎风、背风、近切向及特殊面元标识。
% 本模块不计算压力、热流、驻点、分离或真实遮挡，也不自动选择物理模型。
%
% 调用顺序：
%   mesh = analyzeLocalGeometry(mesh,cfg);
%   freestream = calcFreestream(cfg.flight);
%   mesh = classifyAerothermalRegion(mesh,freestream,cfg);
%
% 输入：mesh 为外法向已修正的最终网格，含 localGeometry；
%       freestream.unitFlowDirection 为相对箭体的空气流动方向；
%       cfg 可省略。cfg.aerothermalRegion.grazingAngleDeg 默认 0.5 度。
% 坐标：鼻尖 +X；零攻角空气沿 -X 流动。不得使用飞行方向代替流动方向。
%
% 输出：新增 mesh.aerothermalRegion，其他字段保持不变。
%   mu = -dot(normal,unitFlowDirection)，正值表示局部迎风；
%   incidenceAngleDeg = asin(mu)，相对切平面的有符号入射角；
%   isWindward / isLeeward / isGrazing 为互斥且完备的三类标识；
%   flowClass 为 windward / leeward / grazing 字符串元胞列；
%   geometryRegion 与端部标识复用 localGeometry，不按攻角重新定义；
%   requiresSpecialTreatment 汇总封盖、极点邻接、折角及网格问题；
%   summary 提供各类面元数量和面积，便于检查，不代表热流占比。
%
% 注意：近切向阈值仅用于分类，原始 mu 保留供后续压力公式使用。
% 凹形体的局部迎风面可能被其他表面遮挡；此处没有射线可见性检测。
% 极点邻接不等于物理尖鼻，迎风也不等于驻点或层流区域。

%% 1. 检查接口：禁止在缺少局部几何或来流方向时猜测输入
if nargin < 3, cfg = struct; end
if nargin < 2 || ~isstruct(freestream) || ~isfield(freestream,'unitFlowDirection')
    error('RocketAeroHeat:Region:Freestream','必须提供 freestream.unitFlowDirection。');
end
if ~isstruct(mesh) || ~all(isfield(mesh,{'V','F','localGeometry'}))
    error('RocketAeroHeat:Region:Geometry','请先对最终网格执行 analyzeLocalGeometry。');
end
flow = freestream.unitFlowDirection;
validateattributes(flow,{'numeric'},{'real','finite','vector','numel',3});
flow = double(flow(:)');
flowNorm = norm(flow);
if ~isfinite(flowNorm) || flowNorm <= 0
    error('RocketAeroHeat:Region:FlowDirection','来流方向必须是非零有限向量。');
end
flow = flow/flowNorm;
threshold = 0.5;
if isfield(cfg,'aerothermalRegion') && isfield(cfg.aerothermalRegion,'grazingAngleDeg')
    threshold = cfg.aerothermalRegion.grazingAngleDeg;
end
validateattributes(threshold,{'numeric'},{'real','finite','scalar','>=',0,'<',90});

%% 2. 检查局部几何与当前面元是否一致，防止编辑网格后使用旧结果
V = mesh.V; F = mesh.F; g = mesh.localGeometry;
validateattributes(V,{'numeric'},{'real','finite','2d','ncols',3,'nonempty'});
validateattributes(F,{'numeric'},{'real','finite','integer','positive','2d','ncols',3,'nonempty','<=',size(V,1)});
V = double(V); F = double(F); nf = size(F,1);
for name = {'normal','center','area','region','isAftCap','isNoseCap', ...
        'isAftTipAdjacent','isNoseTipAdjacent','poorAngle','featureFace','topologyIssue'}
    key = name{1};
    if ~isfield(g,key) || size(g.(key),1) ~= nf
        error('RocketAeroHeat:Region:Geometry','localGeometry.%s 缺失或面元数量不匹配。',key);
    end
end
validateattributes(g.normal,{'numeric'},{'real','finite','size',[nf 3]});
validateattributes(g.center,{'numeric'},{'real','finite','size',[nf 3]});
validateattributes(g.area,{'numeric'},{'real','finite','positive','size',[nf 1]});
for name = {'isAftCap','isNoseCap','isAftTipAdjacent','isNoseTipAdjacent', ...
        'poorAngle','featureFace','topologyIssue'}
    validateattributes(g.(name{1}),{'logical','numeric'},{'binary','size',[nf 1]});
end
p1 = V(F(:,1),:); e1 = V(F(:,2),:)-p1; e2 = V(F(:,3),:)-p1;
raw = cross(e1,e2,2); magnitude = vecnorm(raw,2,2);
if any(magnitude <= 0)
    error('RocketAeroHeat:Region:Geometry','当前网格包含退化面元。');
end
normal = raw./magnitude;
center = p1+(e1+e2)/3;
tol = max(1e-10*max(max(V)-min(V)),64*eps(max(abs(V(:)))));
if any(vecnorm(normal-g.normal,2,2) > 1e-10) || ...
        any(vecnorm(center-g.center,2,2) > tol) || ...
        any(abs(magnitude/2-g.area) > 1e-10*(magnitude/2))
    error('RocketAeroHeat:Region:StaleGeometry','网格已改变，请重新执行局部几何分析。');
end

%% 3. 按法向与空气流动方向分类：三种类别互斥，保留原始入射量
mu = max(-1,min(1,-normal*flow'));
% 微小数值裕量使恰好位于阈值上的面元稳定归入近切向类。
muTolerance = sind(threshold)+32*eps;
windward = mu > muTolerance;
leeward = mu < -muTolerance;
grazing = ~(windward|leeward);
flowClass = repmat({'grazing'},nf,1);
flowClass(windward) = {'windward'};
flowClass(leeward) = {'leeward'};

%% 4. 复用端部和诊断标识：几何类别与来流类别独立，避免覆盖信息
r = struct;
r.version = '1.0'; r.unitFlowDirection = flow;
r.grazingAngleDeg = threshold; r.muTolerance = muTolerance;
r.mu = mu; r.incidenceAngleDeg = asind(mu);
r.isWindward = windward; r.isLeeward = leeward; r.isGrazing = grazing;
r.flowClass = flowClass; r.geometryRegion = g.region;
for name = {'isAftCap','isNoseCap','isAftTipAdjacent','isNoseTipAdjacent', ...
        'poorAngle','featureFace','topologyIssue'}
    r.(name{1}) = logical(g.(name{1}));
end
r.isSidewall = ~(r.isAftCap|r.isNoseCap);
r.requiresSpecialTreatment = r.isAftCap | r.isNoseCap | ...
    r.isAftTipAdjacent | r.isNoseTipAdjacent | r.featureFace | r.poorAngle | r.topologyIssue;
r.classificationBasis = '仅依据局部法向和来流；未计算遮挡、分离、驻点或流态';

%% 5. 汇总分类数量和面积，供主程序与后续可视化读取
r.summary.numberFaces = nf;
r.summary.totalArea = sum(g.area);
for name = {'Windward','Leeward','Grazing'}
    key = name{1}; mask = r.(['is' key]);
    r.summary.([lower(key) 'Faces']) = nnz(mask);
    r.summary.([lower(key) 'Area']) = sum(g.area(mask));
end
r.summary.specialTreatmentFaces = nnz(r.requiresSpecialTreatment);
mesh.aerothermalRegion = r;
verbose = true;
if isfield(cfg,'output') && isfield(cfg.output,'verbose'), verbose = cfg.output.verbose; end
if verbose
    fprintf('\n区域分类完成：迎风 %d，背风 %d，近切向 %d。\n', ...
        nnz(windward),nnz(leeward),nnz(grazing));
    fprintf('近切向角阈值 %.3f 度；特殊处理面元 %d。\n',threshold,nnz(r.requiresSpecialTreatment));
end
end
