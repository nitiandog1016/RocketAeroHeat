function mesh = autoRemeshSurface(rawMesh, cfg)
%% ========================================================================
% autoRemeshSurface.m
%
% Rocket AeroHeat
%
% 高性能版：
% 自适应轴向截面重构半结构化三角表面网格
%
% ========================================================================
%
% 核心网格方法保持不变：
%
%       原始 STL
%           ↓
%       X = constant 截面
%           ↓
%       三角形 / 平面求交
%           ↓
%       Y-Z 闭合轮廓
%           ↓
%       自适应轴向加密
%           ↓
%       周向等弧长重采样
%           ↓
%       Ring-to-Ring zipper triangulation
%           ↓
%       两端封闭
%           ↓
%       初步外法向修正
%
%
% ========================================================================
% 本版本性能优化
%
% 1. Triangle X-range acceleration
%
% 原算法每一个截面都遍历全部 STL 三角形：
%
%       Nsection × Nface
%
% 现在预计算：
%
%       xmin_face
%       xmax_face
%
% 对 x = xs 截面，只检查：
%
%       xmin_face <= xs <= xmax_face
%
%
% 2. Adaptive section cache
%
% 原算法：
%
%       每次加密以后
%       所有旧截面全部重新计算
%
% 现在：
%
%       已经计算过的截面直接从 cache 读取
%
%       只计算新增 midpoint
%
%
% 3. 减少动态数组扩展
%
%       segments
%       faces
%
% 均尽量预分配。
%
%
% ========================================================================
% 坐标约定
%
%       主轴 = X
%       鼻尖 = +X
%       尾部 = -X
%
%
% ========================================================================
% 输出：
%
% mesh.V
% mesh.F
%
% mesh.center
% mesh.normal
% mesh.area
%
% mesh.param.*
%
% mesh.remeshInfo.*
%
% mesh.remeshInfo.timing.*
%
%% ========================================================================



%% ========================================================================
% 1. 总计时
% ========================================================================

timerTotal = tic;



%% ========================================================================
% 2. 输入检查
% ========================================================================

if nargin < 2

    error( ...
        'RocketAeroHeat:Remesh:MissingInput', ...
        'autoRemeshSurface 必须输入 rawMesh 和 cfg。');

end


if ~isstruct(rawMesh) ...
        || ...
        ~isfield(rawMesh,'V') ...
        || ...
        ~isfield(rawMesh,'F')

    error( ...
        'RocketAeroHeat:Remesh:InvalidRawMesh', ...
        'rawMesh 必须包含 V 和 F。');

end


Vraw = ...
    double(rawMesh.V);


Fraw = ...
    double(rawMesh.F);


validateXYZMeshArrays( ...
    Vraw, ...
    Fraw, ...
    'rawMesh');


if any(~isfinite(Vraw(:)))

    error( ...
        'RocketAeroHeat:Remesh:NonFiniteVertex', ...
        'rawMesh.V 中存在 NaN 或 Inf。');

end


if any(~isfinite(Fraw(:))) ...
        || ...
        any(Fraw(:) ~= round(Fraw(:))) ...
        || ...
        any(Fraw(:) < 1) ...
        || ...
        any(Fraw(:) > size(Vraw,1))

    error( ...
        'RocketAeroHeat:Remesh:InvalidFaceIndex', ...
        'rawMesh.F 中存在非法节点编号。');

end



%% ========================================================================
% 3. 参数
% ========================================================================

p = ...
    getParameters(cfg);



%% ========================================================================
% 4. 几何尺度
% ========================================================================

xMin = ...
    min(Vraw(:,1));


xMax = ...
    max(Vraw(:,1));


bodyLength = ...
    xMax-xMin;


if ~(isfinite(bodyLength) ...
        && ...
        bodyLength > 0)

    error( ...
        'RocketAeroHeat:Remesh:ZeroBodyLength', ...
        '几何 X 方向长度必须大于 0。');

end


modelRange = ...
    max(Vraw,[],1) ...
    - ...
    min(Vraw,[],1);


referenceLength = ...
    max(modelRange);


if ~(isfinite(referenceLength) ...
        && ...
        referenceLength > 0)

    referenceLength = ...
        bodyLength;

end



%% ========================================================================
% 5. 建立三角形-截面求交加速数据
%
% 这一部分只计算一次。
%
% 后面所有 x = constant 截面均复用。
%% ========================================================================

timerAcceleration = tic;


accel = ...
    buildIntersectionAccelerator( ...
        Vraw, ...
        Fraw);


timeAcceleration = ...
    toc(timerAcceleration);



%% ========================================================================
% 6. 自适应建立轴向截面
% ========================================================================

fprintf('\n');

fprintf('------------------------------------------------------------\n');
fprintf(' Adaptive Axial Section Remeshing - Performance Version\n');
fprintf('------------------------------------------------------------\n');


timerSections = tic;


[sections,sectionCacheInfo] = ...
    buildAdaptiveSections( ...
        accel, ...
        p);


timeSections = ...
    toc(timerSections);



%% ========================================================================
% 7. 截面数量
% ========================================================================

numberSections = ...
    numel(sections);


if numberSections < 2

    error( ...
        'RocketAeroHeat:Remesh:TooFewSections', ...
        ['有效轴向截面少于 2 个。' newline ...
         '无法建立表面网格。']);

end



%% ========================================================================
% 8. 周向节点数量
% ========================================================================

timerCircumference = tic;


sectionPointCount = ...
    zeros(numberSections,1);


sectionTargetSpacing = ...
    zeros(numberSections,1);


for i = 1:numberSections

    sectionTargetSpacing(i) = ...
        localSectionTargetSpacing( ...
            sections(i).metrics, ...
            p);


    sectionPointCount(i) = ...
        determineCircumferentialPointNumber( ...
            sections(i).metrics, ...
            sectionTargetSpacing(i), ...
            p);

end



%% ========================================================================
% 9. 限制相邻截面周向节点数变化
% ========================================================================

sectionPointCount = ...
    smoothCircumferentialCounts( ...
        sectionPointCount, ...
        p);



%% ========================================================================
% 10. 各截面重采样
% ========================================================================

rings = ...
    cell(numberSections,1);


for i = 1:numberSections

    rings{i} = ...
        resampleClosedLoop( ...
            sections(i).loop, ...
            sectionPointCount(i));


    if size(rings{i},2) ~= 3

        error( ...
            'RocketAeroHeat:Remesh:InvalidRingDimension', ...
            ['第 %d 个截面的重采样结果不是 N×3。' newline ...
             '实际尺寸 = %d×%d'], ...
            i, ...
            size(rings{i},1), ...
            size(rings{i},2));

    end

end


timeCircumference = ...
    toc(timerCircumference);



%% ========================================================================
% 11. 创建所有环节点
% ========================================================================

timerTriangulation = tic;


totalRingNodes = ...
    sum(sectionPointCount);


V = ...
    zeros(totalRingNodes,3);


nodeSectionID = ...
    zeros(totalRingNodes,1);


nodeCircumID = ...
    zeros(totalRingNodes,1);


nodeRingFraction = ...
    zeros(totalRingNodes,1);


nodeIsCap = ...
    false(totalRingNodes,1);


sectionNodeIDs = ...
    cell(numberSections,1);


sectionStartNode = ...
    zeros(numberSections,1);



%% ========================================================================
% 12. 写入环节点
% ========================================================================

currentNode = ...
    0;


for i = 1:numberSections

    ring = ...
        rings{i};


    numberPoint = ...
        size(ring,1);


    nodeIDs = ...
        ( ...
        currentNode+1 : ...
        currentNode+numberPoint ...
        )';


    V(nodeIDs,:) = ...
        ring;


    nodeSectionID(nodeIDs) = ...
        i;


    nodeCircumID(nodeIDs) = ...
        (1:numberPoint)';


    nodeRingFraction(nodeIDs) = ...
        (0:numberPoint-1)' ...
        / ...
        numberPoint;


    nodeIsCap(nodeIDs) = ...
        false;


    sectionNodeIDs{i} = ...
        nodeIDs;


    sectionStartNode(i) = ...
        nodeIDs(1);


    currentNode = ...
        currentNode ...
        + ...
        numberPoint;

end



%% ========================================================================
% 13. 预估侧壁三角面数量
%
% zipper 两个 ring 之间：
%
% 面数量大约：
%
%       NA + NB
%
% 相等节点数时实际是：
%
%       2*N
%% ========================================================================

estimatedSideFaces = ...
    0;


for i = 1:numberSections-1

    estimatedSideFaces = ...
        estimatedSideFaces ...
        + ...
        sectionPointCount(i) ...
        + ...
        sectionPointCount(i+1);

end


estimatedTotalFaces = estimatedSideFaces + sectionPointCount(1) + sectionPointCount(end);
if estimatedTotalFaces > p.maxFaces
    error('RocketAeroHeat:Remesh:FaceLimit', ...
        'Requested mesh needs %d faces, exceeding maxFaces = %d. Increase spacing or the limit.', ...
        estimatedTotalFaces,p.maxFaces);
end
F = ...
    zeros(estimatedSideFaces,3);


faceBandID = ...
    zeros(estimatedSideFaces,1);


currentFace = ...
    0;



%% ========================================================================
% 14. Ring-to-Ring triangulation
% ========================================================================

for i = 1:numberSections-1

    ringA = ...
        sectionNodeIDs{i};


    ringB = ...
        sectionNodeIDs{i+1};


    bandFaces = ...
        stitchTwoRings( ...
            ringA, ...
            ringB);


    numberBandFaces = ...
        size(bandFaces,1);


    targetFaces = ...
        currentFace ...
        + ...
        (1:numberBandFaces);


    F(targetFaces,:) = ...
        bandFaces;


    faceBandID(targetFaces) = ...
        i;


    currentFace = ...
        currentFace ...
        + ...
        numberBandFaces;

end



%% ========================================================================
% 删除预分配后未使用部分
% ========================================================================

F = ...
    F(1:currentFace,:);


faceBandID = ...
    faceBandID(1:currentFace);



%% ========================================================================
% 15. 尾端封闭
% ========================================================================

[V, ...
 aftFaces, ...
 aftCapNodeID, ...
 aftCapType, ...
 nodeSectionID, ...
 nodeCircumID, ...
 nodeRingFraction, ...
 nodeIsCap] = ...
    closeBodyEnd( ...
        V, ...
        rawMesh, ...
        sectionNodeIDs{1}, ...
        'aft', ...
        p, ...
        nodeSectionID, ...
        nodeCircumID, ...
        nodeRingFraction, ...
        nodeIsCap);


F = [
    F;
    aftFaces
];


faceBandID = [
    faceBandID;
    -ones(size(aftFaces,1),1)
];



%% ========================================================================
% 16. 鼻端封闭
% ========================================================================

[V, ...
 noseFaces, ...
 noseCapNodeID, ...
 noseCapType, ...
 nodeSectionID, ...
 nodeCircumID, ...
 nodeRingFraction, ...
 nodeIsCap] = ...
    closeBodyEnd( ...
        V, ...
        rawMesh, ...
        sectionNodeIDs{end}, ...
        'nose', ...
        p, ...
        nodeSectionID, ...
        nodeCircumID, ...
        nodeRingFraction, ...
        nodeIsCap);


F = [
    F;
    noseFaces
];


faceBandID = [
    faceBandID;
    numberSections ...
    * ...
    ones(size(noseFaces,1),1)
];



%% ========================================================================
% 17. 网格维度检查
% ========================================================================

validateXYZMeshArrays( ...
    V, ...
    F, ...
    'generated mesh');


if numel(faceBandID) ~= size(F,1)

    error( ...
        'RocketAeroHeat:Remesh:FaceMetadataMismatch', ...
        'faceBandID 数量与三角面数量不一致。');

end



%% ========================================================================
% 18. 截面中心
% ========================================================================

sectionX = ...
    reshape( ...
        [sections.x], ...
        [], ...
        1);


sectionCenterYZ = ...
    zeros(numberSections,2);


for i = 1:numberSections

    centerYZ = ...
        sections(i).metrics.centerYZ;


    if numel(centerYZ) ~= 2

        error( ...
            'RocketAeroHeat:Remesh:InvalidSectionCenter', ...
            '第 %d 个截面的 centerYZ 必须包含两个元素。', ...
            i);

    end


    sectionCenterYZ(i,:) = ...
        reshape(centerYZ,1,2);

end



%% ========================================================================
% 19. 初步外法向
% ========================================================================

F = ...
    orientFacesOutward( ...
        V, ...
        F, ...
        faceBandID, ...
        sectionX, ...
        sectionCenterYZ, ...
        numberSections);



%% ========================================================================
% 20. Metadata
% ========================================================================

paramWorking = struct;


paramWorking.nodeSectionID = ...
    nodeSectionID;


paramWorking.nodeCircumID = ...
    nodeCircumID;


paramWorking.nodeRingFraction = ...
    nodeRingFraction;


paramWorking.nodeIsCap = ...
    nodeIsCap;


paramWorking.faceBandID = ...
    faceBandID;


paramWorking.aftCapNodeID = ...
    aftCapNodeID;


paramWorking.noseCapNodeID = ...
    noseCapNodeID;



%% ========================================================================
% 21. 清理网格
% ========================================================================

[V, ...
 F, ...
 paramWorking] = ...
    cleanGeneratedMesh( ...
        V, ...
        F, ...
        paramWorking, ...
        referenceLength);


validateXYZMeshArrays( ...
    V, ...
    F, ...
    'cleaned generated mesh');



%% ========================================================================
% 22. 清理后重建 sectionNodeIDs
% ========================================================================

sectionNodeIDs = ...
    cell(numberSections,1);


sectionStartNode = ...
    zeros(numberSections,1);


sectionPointCountFinal = ...
    zeros(numberSections,1);


for i = 1:numberSections

    ids = ...
        find( ...
            paramWorking.nodeSectionID == i);


    sectionNodeIDs{i} = ...
        ids;


    sectionPointCountFinal(i) = ...
        numel(ids);


    if ~isempty(ids)

        sectionStartNode(i) = ...
            ids(1);

    end

end


timeTriangulation = ...
    toc(timerTriangulation);



%% ========================================================================
% 23. Topology / Quality
% ========================================================================

timerQuality = tic;


topology = ...
    checkSurfaceTopology(F);


quality = ...
    calculateTriangleQuality( ...
        V, ...
        F);


timeQuality = ...
    toc(timerQuality);



%% ========================================================================
% 24. 创建输出
% ========================================================================

mesh = struct;


mesh.V = ...
    V;


mesh.F = ...
    F;


mesh = ...
    calculateBasicFaceGeometry(mesh);



%% ========================================================================
% 25. 保存半结构化参数信息
% ========================================================================

mesh.param = struct;


mesh.param.version = ...
    4;


mesh.param.description = ...
    ['Accelerated adaptive axial-section-based ', ...
     'semi-structured triangular surface mesh'];


mesh.param.bodyAxis = ...
    [1 0 0];


mesh.param.noseDirection = ...
    [1 0 0];


mesh.param.nodeSectionID = ...
    paramWorking.nodeSectionID;


mesh.param.nodeCircumID = ...
    paramWorking.nodeCircumID;


mesh.param.nodeRingFraction = ...
    paramWorking.nodeRingFraction;


mesh.param.nodeIsCap = ...
    paramWorking.nodeIsCap;


mesh.param.sectionX = ...
    sectionX;


mesh.param.sectionPointCount = ...
    sectionPointCountFinal;


mesh.param.sectionStartNode = ...
    sectionStartNode;


mesh.param.sectionNodeIDs = ...
    sectionNodeIDs;


mesh.param.sectionTargetSpacing = ...
    sectionTargetSpacing;


mesh.param.sectionPerimeter = ...
    reshape( ...
        arrayfun( ...
            @(s)s.metrics.perimeter, ...
            sections), ...
        [], ...
        1);


mesh.param.sectionAreaYZ = ...
    reshape( ...
        arrayfun( ...
            @(s)s.metrics.areaYZ, ...
            sections), ...
        [], ...
        1);


mesh.param.sectionEquivalentRadius = ...
    reshape( ...
        arrayfun( ...
            @(s)s.metrics.equivalentRadius, ...
            sections), ...
        [], ...
        1);


mesh.param.sectionCenterYZ = ...
    sectionCenterYZ;


mesh.param.faceBandID = ...
    paramWorking.faceBandID;


mesh.param.aftCapNodeID = ...
    paramWorking.aftCapNodeID;


mesh.param.noseCapNodeID = ...
    paramWorking.noseCapNodeID;


mesh.param.aftCapType = ...
    aftCapType;


mesh.param.noseCapType = ...
    noseCapType;


mesh.param.xMinOriginal = ...
    xMin;


mesh.param.xMaxOriginal = ...
    xMax;



%% ========================================================================
% 26. Remesh 信息
% ========================================================================

mesh.remeshInfo = struct;


mesh.remeshInfo.numberSections = ...
    numberSections;

mesh.remeshInfo.minCircumPoints = min(sectionPointCountFinal);
mesh.remeshInfo.maxCircumPoints = max(sectionPointCountFinal);


mesh.remeshInfo.numberVertices = ...
    size(V,1);


mesh.remeshInfo.numberFaces = ...
    size(F,1);


mesh.remeshInfo.boundaryEdges = ...
    topology.boundaryEdges;


mesh.remeshInfo.manifoldEdges = ...
    topology.manifoldEdges;


mesh.remeshInfo.nonManifoldEdges = ...
    topology.nonManifoldEdges;


mesh.remeshInfo.isClosed = ...
    topology.isClosed;


mesh.remeshInfo.meanAspectRatio = ...
    mean(quality.aspectRatio);


mesh.remeshInfo.maxAspectRatio = ...
    max(quality.aspectRatio);


mesh.remeshInfo.meanMinimumAngle = ...
    mean(quality.minAngle);


mesh.remeshInfo.minimumAngle = ...
    min(quality.minAngle);


mesh.remeshInfo.aftCapType = ...
    aftCapType;


mesh.remeshInfo.noseCapType = ...
    noseCapType;



%% ========================================================================
% 27. 性能统计
% ========================================================================

timeTotal = ...
    toc(timerTotal);


mesh.remeshInfo.timing = struct;


mesh.remeshInfo.timing.accelerationSetup = ...
    timeAcceleration;


mesh.remeshInfo.timing.sectionExtraction = ...
    timeSections;


mesh.remeshInfo.timing.circumferentialResampling = ...
    timeCircumference;


mesh.remeshInfo.timing.triangulation = ...
    timeTriangulation;


mesh.remeshInfo.timing.topologyAndQuality = ...
    timeQuality;


mesh.remeshInfo.timing.total = ...
    timeTotal;


mesh.remeshInfo.sectionCache = ...
    sectionCacheInfo;



%% ========================================================================
% 28. 输出
% ========================================================================

fprintf('Number of axial sections : %d\n', ...
    numberSections);


fprintf('Number of vertices       : %d\n', ...
    size(V,1));


fprintf('Number of triangles      : %d\n', ...
    size(F,1));


fprintf('\n');


fprintf('Circumferential points   : %d ~ %d\n', ...
    min(sectionPointCountFinal), ...
    max(sectionPointCountFinal));


fprintf('Mean Aspect Ratio        : %.4f\n', ...
    mean(quality.aspectRatio));


fprintf('Maximum Aspect Ratio     : %.4f\n', ...
    max(quality.aspectRatio));


fprintf('Mean minimum angle       : %.4f deg\n', ...
    mean(quality.minAngle));


fprintf('Minimum triangle angle   : %.4f deg\n', ...
    min(quality.minAngle));


fprintf('\n');


fprintf('Boundary edges           : %d\n', ...
    topology.boundaryEdges);


fprintf('Non-manifold edges       : %d\n', ...
    topology.nonManifoldEdges);


fprintf('Aft end                  : %s\n', ...
    aftCapType);


fprintf('Nose end                 : %s\n', ...
    noseCapType);


fprintf('\n');


fprintf('Performance statistics:\n');


fprintf('  Acceleration setup     : %.4f s\n', ...
    timeAcceleration);


fprintf('  Section extraction     : %.4f s\n', ...
    timeSections);


fprintf('  Ring resampling        : %.4f s\n', ...
    timeCircumference);


fprintf('  Triangulation          : %.4f s\n', ...
    timeTriangulation);


fprintf('  Topology / quality     : %.4f s\n', ...
    timeQuality);


fprintf('  Total remesh time      : %.4f s\n', ...
    timeTotal);


fprintf('\n');


fprintf('Section cache:\n');


fprintf('  Actual section solves  : %d\n', ...
    sectionCacheInfo.numberComputed);


fprintf('  Cached section reuses  : %d\n', ...
    sectionCacheInfo.numberReused);


fprintf('------------------------------------------------------------\n');

end



%% ========================================================================
% LOCAL FUNCTION
% 建立截面求交加速器
%% ========================================================================
function accel = ...
    buildIntersectionAccelerator(V,F)

validateXYZMeshArrays( ...
    V, ...
    F, ...
    'intersection accelerator');


%% ------------------------------------------------------------------------
% 每个三角形三个顶点
%% ------------------------------------------------------------------------

P1 = ...
    V(F(:,1),:);


P2 = ...
    V(F(:,2),:);


P3 = ...
    V(F(:,3),:);



%% ------------------------------------------------------------------------
% 每个三角形 X 范围
%% ------------------------------------------------------------------------

faceXMin = ...
    min( ...
        [ ...
        P1(:,1), ...
        P2(:,1), ...
        P3(:,1) ...
        ], ...
        [], ...
        2);


faceXMax = ...
    max( ...
        [ ...
        P1(:,1), ...
        P2(:,1), ...
        P3(:,1) ...
        ], ...
        [], ...
        2);



%% ------------------------------------------------------------------------
% 保存
%% ------------------------------------------------------------------------

accel = struct;


accel.V = ...
    V;


accel.F = ...
    F;


accel.P1 = ...
    P1;


accel.P2 = ...
    P2;


accel.P3 = ...
    P3;


accel.faceXMin = ...
    faceXMin;


accel.faceXMax = ...
    faceXMax;


accel.numberFaces = ...
    size(F,1);


accel.xMin = ...
    min(V(:,1));


accel.xMax = ...
    max(V(:,1));


accel.modelRange = ...
    max(V,[],1) ...
    - ...
    min(V,[],1);


accel.referenceLength = ...
    max(accel.modelRange);

end



%% ========================================================================
% LOCAL FUNCTION
% 参数
%% ========================================================================
function p = getParameters(cfg)

p.hMin = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','hMin'}, ...
        0.003);


p.hMax = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','hMax'}, ...
        0.030);


p.initialAxialSpacing = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','initialAxialSpacing'}, ...
        0.030);


p.axialSpacingFactor = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','axialSpacingFactor'}, ...
        1.25);


p.maxSectionPerimeterChange = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','maxSectionPerimeterChange'}, ...
        0.12);


p.sectionRefineIterations = ...
    round( ...
        getFieldOrDefault( ...
            cfg, ...
            {'remesh','sectionRefineIterations'}, ...
            8));


p.maxStations = ...
    round( ...
        getFieldOrDefault( ...
            cfg, ...
            {'remesh','maxStations'}, ...
            1200));


p.curvatureAngleTarget = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','curvatureAngleTarget'}, ...
        8.0);


p.curvatureSensitivity = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','curvatureSensitivity'}, ...
        1.0);


p.minCircumPoints = ...
    round( ...
        getFieldOrDefault( ...
            cfg, ...
            {'remesh','minCircumPoints'}, ...
            8));


p.maxCircumPoints = ...
    round( ...
        getFieldOrDefault( ...
            cfg, ...
            {'remesh','maxCircumPoints'}, ...
            160));


p.maxCircumCountRatio = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','maxCircumCountRatio'}, ...
        1.35);


p.maxAspectRatio = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','maxAspectRatio'}, ...
        5.0);


p.minAngle = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','minAngle'}, ...
        20.0);


p.sectionMergeToleranceFraction = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','sectionMergeToleranceFraction'}, ...
        1e-7);


p.planeToleranceFraction = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','planeToleranceFraction'}, ...
        1e-10);


p.minimumEndOffsetFraction = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','minimumEndOffsetFraction'}, ...
        1e-5);


p.flatEndRadiusRatio = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','flatEndRadiusRatio'}, ...
        0.40);


p.extremeVertexToleranceFraction = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','extremeVertexToleranceFraction'}, ...
        1e-6);


p.minFaces = ...
    round( ...
        getFieldOrDefault( ...
            cfg, ...
            {'remesh','minFaces'}, ...
            3000));


p.maxFaces = ...
    round( ...
        getFieldOrDefault( ...
            cfg, ...
            {'remesh','maxFaces'}, ...
            100000));


p.targetFaceNumber = ...
    getFieldOrDefault( ...
        cfg, ...
        {'remesh','targetFaceNumber'}, ...
        []);



%% ------------------------------------------------------------------------
% 合法性
%% ------------------------------------------------------------------------

positiveFields = {'hMin','hMax','initialAxialSpacing','axialSpacingFactor', ...
    'maxSectionPerimeterChange','curvatureAngleTarget','curvatureSensitivity', ...
    'sectionMergeToleranceFraction','planeToleranceFraction', ...
    'minimumEndOffsetFraction','flatEndRadiusRatio','extremeVertexToleranceFraction'};
for k = 1:numel(positiveFields)
    name = positiveFields{k};
    validateattributes(p.(name), {'numeric'}, ...
        {'scalar','real','finite','positive'},mfilename,name);
end
validateattributes(p.maxStations, {'numeric'}, {'scalar','real','finite','integer','>=',3});
validateattributes(p.sectionRefineIterations, {'numeric'}, {'scalar','real','finite','integer','>=',0});
validateattributes(p.minCircumPoints, {'numeric'}, {'scalar','real','finite','integer','>=',3});
validateattributes(p.maxCircumPoints, {'numeric'}, ...
    {'scalar','real','finite','integer','>=',p.minCircumPoints});
validateattributes(p.maxFaces, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(p.maxCircumCountRatio, {'numeric'}, {'scalar','real','finite','>=',1});
validateattributes(p.minimumEndOffsetFraction, {'numeric'}, {'<',0.5});

if p.hMin <= 0 ...
        || ...
        p.hMax <= 0 ...
        || ...
        p.hMin >= p.hMax

    error( ...
        'RocketAeroHeat:Remesh:InvalidSizeRange', ...
        '必须满足 0 < hMin < hMax。');

end


if p.initialAxialSpacing <= 0

    error( ...
        'RocketAeroHeat:Remesh:InvalidAxialSpacing', ...
        'initialAxialSpacing 必须大于 0。');

end


p.minCircumPoints = ...
    max(3,p.minCircumPoints);


p.maxCircumPoints = ...
    max( ...
        p.minCircumPoints, ...
        p.maxCircumPoints);


p.maxStations = ...
    max(3,p.maxStations);


p.sectionRefineIterations = ...
    max(0,p.sectionRefineIterations);

end



%% ========================================================================
% LOCAL FUNCTION
% cfg 字段
%% ========================================================================
function value = ...
    getFieldOrDefault(cfg,path,defaultValue)

value = ...
    defaultValue;


current = ...
    cfg;


for i = 1:numel(path)

    fieldName = ...
        path{i};


    if ~isstruct(current) ...
            || ...
            ~isfield(current,fieldName)

        return;

    end


    current = ...
        current.(fieldName);

end


value = ...
    current;

end



%% ========================================================================
% LOCAL FUNCTION
% 自适应轴向截面 + Cache
%% ========================================================================
function [sections,cacheInfo] = ...
    buildAdaptiveSections(accel,p)

xMin = ...
    accel.xMin;


xMax = ...
    accel.xMax;


bodyLength = ...
    xMax-xMin;


referenceLength = ...
    accel.referenceLength;


if referenceLength <= 0

    referenceLength = ...
        bodyLength;

end



%% ------------------------------------------------------------------------
% 避开精确端点
%% ------------------------------------------------------------------------

endOffset = ...
    max( ...
        p.minimumEndOffsetFraction ...
        * ...
        bodyLength, ...
        100 ...
        * ...
        eps(referenceLength));


xStart = ...
    xMin ...
    + ...
    endOffset;


xEnd = ...
    xMax ...
    - ...
    endOffset;


if xStart >= xEnd

    xStart = ...
        xMin ...
        + ...
        0.01*bodyLength;


    xEnd = ...
        xMax ...
        - ...
        0.01*bodyLength;

end



%% ------------------------------------------------------------------------
% 初始轴向截面
%% ------------------------------------------------------------------------

numberInitial = ...
    ceil( ...
        (xEnd-xStart) ...
        / ...
        p.initialAxialSpacing) ...
    + ...
    1;


numberInitial = ...
    max(3,numberInitial);


numberInitial = ...
    min( ...
        numberInitial, ...
        p.maxStations);


xStation = ...
    linspace( ...
        xStart, ...
        xEnd, ...
        numberInitial)';



%% ========================================================================
% Cache
%% ========================================================================

emptySection = ...
    createEmptySection();


cacheX = ...
    zeros(0,1);


cacheSections = ...
    repmat( ...
        emptySection, ...
        0, ...
        1);


numberComputed = ...
    0;


numberReused = ...
    0;



%% ========================================================================
% 自适应加密
%% ========================================================================

for iteration = 1:p.sectionRefineIterations


    [sectionsNow, ...
     cacheX, ...
     cacheSections, ...
     computedNow, ...
     reusedNow] = ...
        extractSectionsCached( ...
            accel, ...
            xStation, ...
            p, ...
            cacheX, ...
            cacheSections);


    numberComputed = ...
        numberComputed ...
        + ...
        computedNow;


    numberReused = ...
        numberReused ...
        + ...
        reusedNow;



    %% --------------------------------------------------------------------
    % 删除无效截面
    %% --------------------------------------------------------------------

    valid = ...
        reshape( ...
            [sectionsNow.valid], ...
            [], ...
            1);


    xStation = ...
        xStation(valid);


    sectionsNow = ...
        sectionsNow(valid);


    if numel(xStation) < 2

        error( ...
            'RocketAeroHeat:Remesh:SectionExtractionFailed', ...
            '有效轴向截面不足 2 个。');

    end


    numberStation = ...
        numel(xStation);


    if numberStation >= p.maxStations

        break;

    end



    %% --------------------------------------------------------------------
    % 判断哪些区间需要加密
    %% --------------------------------------------------------------------

    addMidpoint = ...
        false(numberStation-1,1);


    refineScore = ...
        zeros(numberStation-1,1);


    for i = 1:numberStation-1

        dx = ...
            xStation(i+1) ...
            - ...
            xStation(i);


        h1 = ...
            localSectionTargetSpacing( ...
                sectionsNow(i).metrics, ...
                p);


        h2 = ...
            localSectionTargetSpacing( ...
                sectionsNow(i+1).metrics, ...
                p);


        hAxial = ...
            p.axialSpacingFactor ...
            * ...
            min(h1,h2);


        perimeter1 = ...
            sectionsNow(i).metrics.perimeter;


        perimeter2 = ...
            sectionsNow(i+1).metrics.perimeter;


        relativePerimeterChange = ...
            abs( ...
                perimeter2 ...
                - ...
                perimeter1) ...
            / ...
            max( ...
                0.5 ...
                * ...
                (perimeter1+perimeter2), ...
                eps);


        spacingScore = ...
            dx ...
            / ...
            max(hAxial,eps);


        perimeterScore = ...
            relativePerimeterChange ...
            / ...
            max( ...
                p.maxSectionPerimeterChange, ...
                eps);


        refineScore(i) = ...
            max( ...
                spacingScore, ...
                perimeterScore);


        addMidpoint(i) = ...
            spacingScore > 1 ...
            || ...
            perimeterScore > 1;

    end



    %% --------------------------------------------------------------------
    % 不需要继续加密
    %% --------------------------------------------------------------------

    if ~any(addMidpoint)

        break;

    end



    %% --------------------------------------------------------------------
    % 新 midpoint
    %% --------------------------------------------------------------------

    midpoint = ...
        0.5 ...
        * ...
        ( ...
        xStation(1:end-1) ...
        + ...
        xStation(2:end));


    selectedMidpoint = ...
        midpoint(addMidpoint);



    %% --------------------------------------------------------------------
    % maxStations
    %% --------------------------------------------------------------------

    availableNumber = ...
        p.maxStations ...
        - ...
        numel(xStation);


    if availableNumber <= 0

        break;

    end


    if numel(selectedMidpoint) > availableNumber

        candidateScore = ...
            refineScore(addMidpoint);


        [~,order] = ...
            sort( ...
                candidateScore, ...
                'descend');


        selectedMidpoint = ...
            selectedMidpoint( ...
                order(1:availableNumber));

    end



    %% --------------------------------------------------------------------
    % 新截面集合
    %% --------------------------------------------------------------------

    newX = [
        xStation;
        selectedMidpoint(:)
    ];


    newX = ...
        unique( ...
            sort(newX));


    if numel(newX) == numel(xStation)

        break;

    end


    xStation = ...
        newX;

end



%% ========================================================================
% 最终截面
%
% 这里绝大多数都会直接从 cache 读取。
%% ========================================================================

[sections, ...
 cacheX, ...
 cacheSections, ...
 computedNow, ...
 reusedNow] = ...
    extractSectionsCached( ...
        accel, ...
        xStation, ...
        p, ...
        cacheX, ...
        cacheSections); %#ok<ASGLU>


numberComputed = ...
    numberComputed ...
    + ...
    computedNow;


numberReused = ...
    numberReused ...
    + ...
    reusedNow;



%% ------------------------------------------------------------------------
% 删除无效截面
%% ------------------------------------------------------------------------

valid = ...
    reshape( ...
        [sections.valid], ...
        [], ...
        1);


sections = ...
    sections(valid);


if isempty(sections)

    error( ...
        'RocketAeroHeat:Remesh:NoFinalSections', ...
        '最终没有得到有效截面。');

end



%% ------------------------------------------------------------------------
% 按 X 排序
%% ------------------------------------------------------------------------

[~,order] = ...
    sort([sections.x]);


sections = ...
    sections(order);



%% ------------------------------------------------------------------------
% Cache 信息
%% ------------------------------------------------------------------------

cacheInfo = struct;


cacheInfo.numberComputed = ...
    numberComputed;


cacheInfo.numberReused = ...
    numberReused;


cacheInfo.numberUniqueCached = ...
    numel(cacheX);

end



%% ========================================================================
% LOCAL FUNCTION
% Cached section extraction
%% ========================================================================
function [sections, ...
          cacheX, ...
          cacheSections, ...
          numberComputed, ...
          numberReused] = ...
    extractSectionsCached( ...
        accel, ...
        xStation, ...
        p, ...
        cacheX, ...
        cacheSections)

xStation = ...
    xStation(:);


numberStation = ...
    numel(xStation);


emptySection = ...
    createEmptySection();


sections = ...
    repmat( ...
        emptySection, ...
        numberStation, ...
        1);


numberComputed = ...
    0;


numberReused = ...
    0;



%% ========================================================================
% 精确匹配已存在 x
%
% 旧截面的 x 值直接继承，不会重新计算，
% 因此这里 exact match 是可靠的。
%% ========================================================================

if isempty(cacheX)

    cached = ...
        false(numberStation,1);


    location = ...
        zeros(numberStation,1);

else

    [cached,location] = ...
        ismember( ...
            xStation, ...
            cacheX);

end



%% ========================================================================
% 读取 cache
%% ========================================================================

cachedIndex = ...
    find(cached);


for k = 1:numel(cachedIndex)

    i = ...
        cachedIndex(k);


    sections(i) = ...
        cacheSections(location(i));

end


numberReused = ...
    numel(cachedIndex);



%% ========================================================================
% 只计算新截面
%% ========================================================================

missingIndex = ...
    find(~cached);


for k = 1:numel(missingIndex)

    i = ...
        missingIndex(k);


    x = ...
        xStation(i);


    section = ...
        extractOneSection( ...
            accel, ...
            x, ...
            p);


    sections(i) = ...
        section;


    cacheX(end+1,1) = ...
        x; %#ok<AGROW>


    cacheSections(end+1,1) = ...
        section; %#ok<AGROW>


    numberComputed = ...
        numberComputed ...
        + ...
        1;

end

end



%% ========================================================================
% LOCAL FUNCTION
% Empty section
%% ========================================================================
function section = ...
    createEmptySection()

section = struct( ...
    'x',[], ...
    'loop',zeros(0,3), ...
    'metrics',[], ...
    'valid',false);

end



%% ========================================================================
% LOCAL FUNCTION
% 计算单个截面
%% ========================================================================
function section = ...
    extractOneSection(accel,x,p)

section = ...
    createEmptySection();


section.x = ...
    x;


[loop,success] = ...
    extractSectionLoop( ...
        accel, ...
        x, ...
        p);


section.valid = ...
    success;

if ~success
    error('RocketAeroHeat:Remesh:InvalidSection', ...
        'Cannot reconstruct a closed section at x = %.12g. Check the STL or section tolerances.',x);
end


if success

    section.loop = ...
        loop;


    section.metrics = ...
        calculateLoopMetrics(loop);

end

end



%% ========================================================================
% LOCAL FUNCTION
% 提取 X = constant 截面
%
% 性能关键：
%
% 不再遍历全部 STL faces。
%% ========================================================================
function [loop,success] = ...
    extractSectionLoop(accel,xPlane,p)

referenceLength = ...
    accel.referenceLength;


planeTolerance = ...
    max( ...
        p.planeToleranceFraction ...
        * ...
        referenceLength, ...
        100 ...
        * ...
        eps(referenceLength));


mergeTolerance = ...
    max( ...
        p.sectionMergeToleranceFraction ...
        * ...
        referenceLength, ...
        10 ...
        * ...
        planeTolerance);



%% ========================================================================
% X-range 快速筛选
%
% 三角形若满足：
%
%       xmax < xPlane
%
% 或：
%
%       xmin > xPlane
%
% 必然不会与截面相交。
%% ========================================================================

candidateMask = ...
    accel.faceXMin ...
    <= ...
    xPlane+planeTolerance ...
    & ...
    accel.faceXMax ...
    >= ...
    xPlane-planeTolerance;


candidateFaces = ...
    find(candidateMask);


numberCandidates = ...
    numel(candidateFaces);


if numberCandidates < 3

    loop = ...
        zeros(0,3);


    success = ...
        false;


    return;

end



%% ========================================================================
% 最多一个 triangle 对当前 plane 贡献一条 segment
%
% 因此预分配：
%
%       numberCandidates × 6
%% ========================================================================

segments = ...
    zeros(numberCandidates,6);


segmentCount = ...
    0;



%% ========================================================================
% 仅遍历候选三角形
%% ========================================================================

for k = 1:numberCandidates

    face = ...
        candidateFaces(k);


    P1 = ...
        accel.P1(face,:);


    P2 = ...
        accel.P2(face,:);


    P3 = ...
        accel.P3(face,:);


    [point1,point2,hit] = ...
        trianglePlaneIntersectionFast( ...
            P1, ...
            P2, ...
            P3, ...
            xPlane, ...
            planeTolerance);


    if hit

        if norm(point2-point1) ...
                > ...
                mergeTolerance*0.1


            segmentCount = ...
                segmentCount ...
                + ...
                1;


            segments(segmentCount,:) = ...
                [point1 point2];

        end

    end

end



%% ========================================================================
% 截掉未使用预分配区域
%% ========================================================================

segments = ...
    segments(1:segmentCount,:);


if segmentCount < 3

    loop = ...
        zeros(0,3);


    success = ...
        false;


    return;

end



%% ========================================================================
% 恢复闭环
%% ========================================================================

[loop,success] = ...
    orderLoopFromEdges( ...
        segments, ...
        mergeTolerance);


if ~success

    loop = ...
        zeros(0,3);


    return;

end


if size(loop,2) ~= 3

    error( ...
        'RocketAeroHeat:Remesh:LoopDimension', ...
        '截面闭环必须为 N×3。');

end



%% ========================================================================
% 精确锁定 X
%% ========================================================================

loop(:,1) = ...
    xPlane;



%% ========================================================================
% +Y 为起点
%% ========================================================================

[~,startIndex] = ...
    max(loop(:,2));


loop = [
    loop(startIndex:end,:);
    loop(1:startIndex-1,:)
];



%% ========================================================================
% 统一 Y-Z 方向
%% ========================================================================

if polygonSignedAreaYZ(loop) < 0

    if size(loop,1) >= 2

        loop = [
            loop(1,:);
            flipud(loop(2:end,:))
        ];

    end

end


success = ...
    size(loop,1) >= 3;

end



%% ========================================================================
% LOCAL FUNCTION
% 高性能 triangle-plane intersection
%% ========================================================================
function [point1,point2,hit] = ...
    trianglePlaneIntersectionFast( ...
        P1, ...
        P2, ...
        P3, ...
        xPlane, ...
        tolerance)

triangle = [
    P1;
    P2;
    P3
];


d = ...
    triangle(:,1) ...
    - ...
    xPlane;

% Coplanar triangles do not define a unique intersection segment. Adjacent
% non-coplanar surface triangles supply the section boundary.
if all(abs(d) <= tolerance)
    point1 = zeros(1,3);
    point2 = zeros(1,3);
    hit = false;
    return;
end



%% ------------------------------------------------------------------------
% 一个三角形最多需要暂存 3 个唯一交点
%% ------------------------------------------------------------------------

points = ...
    zeros(3,3);


pointCount = ...
    0;


edgeList = [
    1 2;
    2 3;
    3 1
];


for e = 1:3

    i = ...
        edgeList(e,1);


    j = ...
        edgeList(e,2);


    Pi = ...
        triangle(i,:);


    Pj = ...
        triangle(j,:);


    di = ...
        d(i);


    dj = ...
        d(j);



    %% --------------------------------------------------------------------
    % 顶点 i 在平面上
    %% --------------------------------------------------------------------

    if abs(di) <= tolerance

        [points,pointCount] = ...
            appendSmallUniquePoint( ...
                points, ...
                pointCount, ...
                Pi, ...
                tolerance);

    end



    %% --------------------------------------------------------------------
    % 顶点 j 在平面上
    %% --------------------------------------------------------------------

    if abs(dj) <= tolerance

        [points,pointCount] = ...
            appendSmallUniquePoint( ...
                points, ...
                pointCount, ...
                Pj, ...
                tolerance);

    end



    %% --------------------------------------------------------------------
    % 边跨越平面
    %% --------------------------------------------------------------------

    if di*dj < 0

        t = ...
            di ...
            / ...
            (di-dj);


        P = ...
            Pi ...
            + ...
            t ...
            * ...
            (Pj-Pi);


        P(1) = ...
            xPlane;


        [points,pointCount] = ...
            appendSmallUniquePoint( ...
                points, ...
                pointCount, ...
                P, ...
                tolerance);

    end

end



%% ========================================================================
% 少于两个交点
%% ========================================================================

if pointCount < 2

    point1 = ...
        zeros(1,3);


    point2 = ...
        zeros(1,3);


    hit = ...
        false;


    return;

end


points = ...
    points(1:pointCount,:);



%% ========================================================================
% 恰好两个点
%% ========================================================================

if pointCount == 2

    point1 = ...
        points(1,:);


    point2 = ...
        points(2,:);


    hit = ...
        norm(point2-point1) ...
        > ...
        tolerance*0.1;


    return;

end



%% ========================================================================
% 三个点：
% 取距离最远的一对
%% ========================================================================

d12 = ...
    sum( ...
        (points(1,:)-points(2,:)).^2);


d23 = ...
    sum( ...
        (points(2,:)-points(3,:)).^2);


d31 = ...
    sum( ...
        (points(3,:)-points(1,:)).^2);


[maximumDistanceSquared,index] = ...
    max( ...
        [d12 d23 d31]);


switch index

    case 1

        point1 = ...
            points(1,:);


        point2 = ...
            points(2,:);


    case 2

        point1 = ...
            points(2,:);


        point2 = ...
            points(3,:);


    otherwise

        point1 = ...
            points(3,:);


        point2 = ...
            points(1,:);

end


hit = ...
    sqrt(maximumDistanceSquared) ...
    > ...
    tolerance*0.1;

end



%% ========================================================================
% LOCAL FUNCTION
% 三角形内部小数组唯一点
%% ========================================================================
function [points,pointCount] = ...
    appendSmallUniquePoint( ...
        points, ...
        pointCount, ...
        P, ...
        tolerance)

P = ...
    reshape(P,1,3);


if pointCount == 0

    pointCount = ...
        1;


    points(1,:) = ...
        P;


    return;

end


existing = ...
    points(1:pointCount,:);


distanceSquared = ...
    sum( ...
        (existing-P).^2, ...
        2);


if all( ...
        distanceSquared ...
        > ...
        tolerance^2)


    if pointCount < size(points,1)

        pointCount = ...
            pointCount ...
            + ...
            1;


        points(pointCount,:) = ...
            P;

    end

end

end



%% ========================================================================
% LOCAL FUNCTION
% 截面线段恢复闭合轮廓
%% ========================================================================
function [loop,success] = ...
    orderLoopFromEdges(segments,tolerance)

success = ...
    false;


loop = ...
    zeros(0,3);


if size(segments,2) ~= 6

    error( ...
        'RocketAeroHeat:Remesh:SegmentDimension', ...
        'segments 必须为 N×6。');

end


numberSegment = ...
    size(segments,1);


allPoints = [
    segments(:,1:3);
    segments(:,4:6)
];


[uniquePoints,pointMap] = ...
    mergePointsGreedy( ...
        allPoints, ...
        tolerance);


if size(uniquePoints,1) < 3

    return;

end



%% ========================================================================
% Edge connectivity
%% ========================================================================

edgeNode1 = ...
    pointMap(1:numberSegment);


edgeNode2 = ...
    pointMap(numberSegment+1:end);


edges = ...
    [edgeNode1 edgeNode2];


edges = ...
    edges( ...
        edges(:,1) ~= edges(:,2), ...
        :);


if isempty(edges)

    return;

end



%% ========================================================================
% 删除重复无向边
%% ========================================================================

edgeSorted = ...
    sort(edges,2);


[~,keepEdge] = ...
    unique( ...
        edgeSorted, ...
        'rows', ...
        'stable');


edges = ...
    edges(keepEdge,:);


numberEdge = ...
    size(edges,1);


numberNode = ...
    size(uniquePoints,1);



%% ========================================================================
% Degree 计数
%% ========================================================================

degree = ...
    accumarray( ...
        edges(:), ...
        1, ...
        [numberNode 1]);

if any(degree(degree > 0) ~= 2)
    error('RocketAeroHeat:Remesh:OpenSection', ...
        'Section edges do not form a simple closed loop. Repair the STL or adjust the merge tolerance.');
end



%% ========================================================================
% 邻接表
%% ========================================================================

adjacency = ...
    cell(numberNode,1);


for e = 1:numberEdge

    a = ...
        edges(e,1);


    b = ...
        edges(e,2);


    adjacency{a}(end+1) = ...
        b;


    adjacency{b}(end+1) = ...
        a;

end



%% ========================================================================
% 连通分量
%% ========================================================================

visited = ...
    false(numberNode,1);


components = ...
    cell(0,1);


for startNode = 1:numberNode

    if visited(startNode) ...
            || ...
            degree(startNode) == 0

        continue;

    end


    queue = ...
        zeros(numberNode,1);


    queueStart = ...
        1;


    queueEnd = ...
        1;


    queue(1) = ...
        startNode;


    visited(startNode) = ...
        true;


    component = ...
        zeros(numberNode,1);


    componentCount = ...
        0;


    while queueStart <= queueEnd

        node = ...
            queue(queueStart);


        queueStart = ...
            queueStart ...
            + ...
            1;


        componentCount = ...
            componentCount ...
            + ...
            1;


        component(componentCount) = ...
            node;


        neighbours = ...
            adjacency{node};


        for k = 1:numel(neighbours)

            nextNode = ...
                neighbours(k);


            if ~visited(nextNode)

                visited(nextNode) = ...
                    true;


                queueEnd = ...
                    queueEnd ...
                    + ...
                    1;


                queue(queueEnd) = ...
                    nextNode;

            end

        end

    end


    if componentCount >= 3

        components{end+1,1} = ...
            component(1:componentCount); %#ok<AGROW>

    end

end


if isempty(components)

    return;

end



%% ========================================================================
% 每个 component 恢复轮廓
%% ========================================================================

if numel(components) ~= 1
    error('RocketAeroHeat:Remesh:MultipleSectionLoops', ...
        'Found %d section loops. This remesher supports a single outer loop; no components were discarded.', ...
        numel(components));
end

candidateLoops = ...
    cell(numel(components),1);


candidatePerimeter = ...
    zeros(numel(components),1);


candidateValid = ...
    false(numel(components),1);


for c = 1:numel(components)

    component = ...
        components{c};


    componentDegree = ...
        degree(component);


    componentPoints = ...
        uniquePoints(component,:);


    componentLoop = ...
        zeros(0,3);



    %% --------------------------------------------------------------------
    % 理想 degree=2 闭环
    %% --------------------------------------------------------------------

    if all(componentDegree == 2)

        [~,localStart] = ...
            max(componentPoints(:,2));


        startNode = ...
            component(localStart);


        currentNode = ...
            startNode;


        previousNode = ...
            0;


        orderedNode = ...
            zeros(numel(component),1);


        orderedCount = ...
            0;


        closed = ...
            false;


        for step = 1:numel(component)+1

            if currentNode == startNode ...
                    && ...
                    step > 1

                closed = ...
                    true;


                break;

            end


            orderedCount = ...
                orderedCount ...
                + ...
                1;


            orderedNode(orderedCount) = ...
                currentNode;


            neighbours = ...
                adjacency{currentNode};


            if previousNode == 0

                nextNode = ...
                    neighbours(1);

            else

                candidate = ...
                    neighbours( ...
                        neighbours ~= previousNode);


                if isempty(candidate)

                    break;

                end


                nextNode = ...
                    candidate(1);

            end


            previousNode = ...
                currentNode;


            currentNode = ...
                nextNode;

        end


        if closed

            orderedNode = ...
                orderedNode(1:orderedCount);


            if numel(orderedNode) >= 3

                componentLoop = ...
                    uniquePoints(orderedNode,:);

            end

        end

    end



    %% --------------------------------------------------------------------
    % Fallback：极角排序
    %% --------------------------------------------------------------------

    if size(componentLoop,1) < 3

        centerYZ = ...
            mean( ...
                componentPoints(:,2:3), ...
                1);


        angle = ...
            atan2( ...
                componentPoints(:,3)-centerYZ(2), ...
                componentPoints(:,2)-centerYZ(1));


        [~,order] = ...
            sort(angle);


        componentLoop = ...
            componentPoints(order,:);

    end


    if size(componentLoop,1) < 3

        continue;

    end



    %% --------------------------------------------------------------------
    % 周长
    %% --------------------------------------------------------------------

    nextPoint = ...
        componentLoop([2:end 1],:);


    perimeter = ...
        sum( ...
            sqrt( ...
                sum( ...
                    (nextPoint-componentLoop).^2, ...
                    2)));


    if perimeter > 0

        candidateLoops{c} = ...
            componentLoop;


        candidatePerimeter(c) = ...
            perimeter;


        candidateValid(c) = ...
            true;

    end

end



%% ========================================================================
% 取最大主闭环
%% ========================================================================

if ~any(candidateValid)

    return;

end


candidatePerimeter(~candidateValid) = ...
    -inf;


[~,largest] = ...
    max(candidatePerimeter);


loop = ...
    candidateLoops{largest};


success = ...
    size(loop,1) >= 3;

end



%% ========================================================================
% LOCAL FUNCTION
% Greedy merge
%% ========================================================================
function [uniquePoints,map] = ...
    mergePointsGreedy(points,tolerance)

if size(points,2) ~= 3

    error( ...
        'RocketAeroHeat:Remesh:MergeInputDimension', ...
        'mergePointsGreedy 要求 N×3。');

end


numberPoint = ...
    size(points,1);


uniquePoints = ...
    zeros(numberPoint,3);


map = ...
    zeros(numberPoint,1);


uniqueCount = ...
    0;


for i = 1:numberPoint

    P = ...
        points(i,:);


    if uniqueCount == 0

        uniqueCount = ...
            1;


        uniquePoints(1,:) = ...
            P;


        map(i) = ...
            1;


        continue;

    end


    existing = ...
        uniquePoints(1:uniqueCount,:);


    distanceSquared = ...
        sum( ...
            (existing-P).^2, ...
            2);


    [minimumDistanceSquared,index] = ...
        min(distanceSquared);


    if minimumDistanceSquared ...
            <= ...
            tolerance^2

        map(i) = ...
            index;


    else

        uniqueCount = ...
            uniqueCount ...
            + ...
            1;


        uniquePoints(uniqueCount,:) = ...
            P;


        map(i) = ...
            uniqueCount;

    end

end


uniquePoints = ...
    uniquePoints(1:uniqueCount,:);

end



%% ========================================================================
% LOCAL FUNCTION
% 截面指标
%% ========================================================================
function metrics = ...
    calculateLoopMetrics(loop)

if size(loop,2) ~= 3 ...
        || ...
        size(loop,1) < 3

    error( ...
        'RocketAeroHeat:Remesh:InvalidMetricLoop', ...
        'calculateLoopMetrics 要求至少 3×3。');

end


nextLoop = ...
    loop([2:end 1],:);


segmentLength = ...
    sqrt( ...
        sum( ...
            (nextLoop-loop).^2, ...
            2));


perimeter = ...
    sum(segmentLength);


signedArea = ...
    polygonSignedAreaYZ(loop);


areaYZ = ...
    abs(signedArea);


centerYZ = ...
    mean(loop(:,2:3),1);


if areaYZ > eps

    equivalentRadius = ...
        sqrt(areaYZ/pi);

else

    equivalentRadius = ...
        perimeter/(2*pi);

end


radialDistance = ...
    sqrt( ...
        (loop(:,2)-centerYZ(1)).^2 ...
        + ...
        (loop(:,3)-centerYZ(2)).^2);


metrics = struct;


metrics.perimeter = ...
    perimeter;


metrics.signedAreaYZ = ...
    signedArea;


metrics.areaYZ = ...
    areaYZ;


metrics.centerYZ = ...
    reshape(centerYZ,1,2);


metrics.equivalentRadius = ...
    equivalentRadius;


metrics.minimumRadius = ...
    min(radialDistance);


metrics.maximumRadius = ...
    max(radialDistance);


metrics.meanRadius = ...
    mean(radialDistance);

end



%% ========================================================================
% LOCAL FUNCTION
% Y-Z 有向面积
%% ========================================================================
function area = ...
    polygonSignedAreaYZ(loop)

y = ...
    loop(:,2);


z = ...
    loop(:,3);


yNext = ...
    y([2:end 1]);


zNext = ...
    z([2:end 1]);


area = ...
    0.5 ...
    * ...
    sum( ...
        y.*zNext ...
        - ...
        yNext.*z);

end



%% ========================================================================
% LOCAL FUNCTION
% 局部周向尺寸
%% ========================================================================
function h = ...
    localSectionTargetSpacing(metrics,p)

radius = ...
    metrics.equivalentRadius;


thetaTarget = ...
    deg2rad( ...
        p.curvatureAngleTarget);


hCurvature = ...
    radius ...
    * ...
    thetaTarget ...
    / ...
    p.curvatureSensitivity;


if ~isfinite(hCurvature) ...
        || ...
        hCurvature <= 0

    hCurvature = ...
        p.hMax;

end


h = ...
    min( ...
        p.hMax, ...
        max( ...
            p.hMin, ...
            hCurvature));

end



%% ========================================================================
% LOCAL FUNCTION
% 周向节点数量
%% ========================================================================
function numberPoint = ...
    determineCircumferentialPointNumber( ...
        metrics, ...
        targetSpacing, ...
        p)

numberPoint = ...
    ceil( ...
        metrics.perimeter ...
        / ...
        max(targetSpacing,eps));


numberPoint = ...
    max( ...
        p.minCircumPoints, ...
        numberPoint);


numberPoint = ...
    min( ...
        p.maxCircumPoints, ...
        numberPoint);


numberPoint = ...
    max( ...
        3, ...
        round(numberPoint));

end



%% ========================================================================
% LOCAL FUNCTION
% 平滑周向节点数量
%% ========================================================================
function counts = ...
    smoothCircumferentialCounts(counts,p)

counts = ...
    round(counts(:));


ratio = ...
    max( ...
        p.maxCircumCountRatio, ...
        1.01);


for iteration = 1:8

    previous = ...
        counts;


    for i = 2:numel(counts)

        maximumAllowed = ...
            ceil( ...
                counts(i-1)*ratio);


        minimumAllowed = ...
            floor( ...
                counts(i-1)/ratio);


        counts(i) = ...
            min( ...
                counts(i), ...
                maximumAllowed);


        counts(i) = ...
            max( ...
                counts(i), ...
                max( ...
                    p.minCircumPoints, ...
                    minimumAllowed));

    end


    for i = numel(counts)-1:-1:1

        maximumAllowed = ...
            ceil( ...
                counts(i+1)*ratio);


        minimumAllowed = ...
            floor( ...
                counts(i+1)/ratio);


        counts(i) = ...
            min( ...
                counts(i), ...
                maximumAllowed);


        counts(i) = ...
            max( ...
                counts(i), ...
                max( ...
                    p.minCircumPoints, ...
                    minimumAllowed));

    end


    counts = ...
        min( ...
            p.maxCircumPoints, ...
            max( ...
                p.minCircumPoints, ...
                counts));


    if isequal(counts,previous)

        break;

    end

end

end



%% ========================================================================
% LOCAL FUNCTION
% 闭合截面等弧长采样
%% ========================================================================
function points = ...
    resampleClosedLoop(loop,numberPoint)

if size(loop,2) ~= 3 ...
        || ...
        size(loop,1) < 3

    error( ...
        'RocketAeroHeat:Remesh:InvalidResampleLoop', ...
        'resampleClosedLoop 要求至少 3×3。');

end


numberPoint = ...
    max( ...
        3, ...
        round(numberPoint));



%% ------------------------------------------------------------------------
% 删除相邻重复点
%% ------------------------------------------------------------------------

nextLoop = ...
    loop([2:end 1],:);


distance = ...
    sqrt( ...
        sum( ...
            (nextLoop-loop).^2, ...
            2));


distanceTolerance = ...
    max( ...
        max(distance)*1e-12, ...
        eps);


keep = ...
    distance ...
    > ...
    distanceTolerance;


loop = ...
    loop(keep,:);


if size(loop,1) < 3

    error( ...
        'RocketAeroHeat:Remesh:TooFewLoopPoints', ...
        '清理重复点后截面不足三个点。');

end



%% ------------------------------------------------------------------------
% CCW
%% ------------------------------------------------------------------------

if polygonSignedAreaYZ(loop) < 0

    loop = ...
        flipud(loop);

end



%% ------------------------------------------------------------------------
% +Y 起点
%% ------------------------------------------------------------------------

[~,startIndex] = ...
    max(loop(:,2));


loop = [
    loop(startIndex:end,:);
    loop(1:startIndex-1,:)
];



%% ------------------------------------------------------------------------
% 闭合
%% ------------------------------------------------------------------------

closedLoop = [
    loop;
    loop(1,:)
];


segment = ...
    diff( ...
        closedLoop, ...
        1, ...
        1);


segmentLength = ...
    sqrt( ...
        sum( ...
            segment.^2, ...
            2));


cumulativeS = [
    0;
    cumsum(segmentLength)
];


totalLength = ...
    cumulativeS(end);


if ~(isfinite(totalLength) ...
        && ...
        totalLength > 0)

    error( ...
        'RocketAeroHeat:Remesh:ZeroPerimeter', ...
        '截面周长为 0 或无效。');

end



%% ------------------------------------------------------------------------
% 等弧长
%% ------------------------------------------------------------------------

sampleS = ...
    linspace( ...
        0, ...
        totalLength, ...
        numberPoint+1)';


sampleS(end) = [];


points = ...
    zeros(numberPoint,3);


for coordinate = 1:3

    points(:,coordinate) = ...
        interp1( ...
            cumulativeS, ...
            closedLoop(:,coordinate), ...
            sampleS, ...
            'linear');

end


points(:,1) = ...
    mean(loop(:,1));


if any(~isfinite(points(:)))

    error( ...
        'RocketAeroHeat:Remesh:NonFiniteResample', ...
        '截面重采样产生 NaN 或 Inf。');

end


if polygonSignedAreaYZ(points) < 0

    points = [
        points(1,:);
        flipud(points(2:end,:))
    ];

end

end



%% ========================================================================
% LOCAL FUNCTION
% Ring zipper
%% ========================================================================
function faces = ...
    stitchTwoRings(ringA,ringB)

ringA = ...
    ringA(:);


ringB = ...
    ringB(:);


nA = ...
    numel(ringA);


nB = ...
    numel(ringB);


if nA < 3 ...
        || ...
        nB < 3

    error( ...
        'RocketAeroHeat:Remesh:TooFewRingNodes', ...
        '相邻两个 ring 都至少需要 3 个节点。');

end



%% ------------------------------------------------------------------------
% 最大三角面数不会超过 nA+nB
%% ------------------------------------------------------------------------

faces = ...
    zeros(nA+nB,3);


faceCount = ...
    0;


stepA = ...
    0;


stepB = ...
    0;


tolerance = ...
    1e-12;


while stepA < nA ...
        || ...
        stepB < nB


    currentA = ...
        ringA( ...
            mod(stepA,nA)+1);


    currentB = ...
        ringB( ...
            mod(stepB,nB)+1);


    if stepA < nA

        nextTA = ...
            (stepA+1)/nA;

    else

        nextTA = ...
            inf;

    end


    if stepB < nB

        nextTB = ...
            (stepB+1)/nB;

    else

        nextTB = ...
            inf;

    end



    %% --------------------------------------------------------------------
    % A / B 同时推进
    %% --------------------------------------------------------------------

    if isfinite(nextTA) ...
            && ...
            isfinite(nextTB) ...
            && ...
            abs(nextTA-nextTB) <= tolerance


        nextA = ...
            ringA( ...
                mod(stepA+1,nA)+1);


        nextB = ...
            ringB( ...
                mod(stepB+1,nB)+1);


        faceCount = ...
            faceCount+1;


        faces(faceCount,:) = ...
            [ ...
            currentA, ...
            nextA, ...
            currentB ...
            ];


        faceCount = ...
            faceCount+1;


        faces(faceCount,:) = ...
            [ ...
            nextA, ...
            nextB, ...
            currentB ...
            ];


        stepA = ...
            stepA+1;


        stepB = ...
            stepB+1;



    %% --------------------------------------------------------------------
    % A 推进
    %% --------------------------------------------------------------------

    elseif nextTA < nextTB


        nextA = ...
            ringA( ...
                mod(stepA+1,nA)+1);


        faceCount = ...
            faceCount+1;


        faces(faceCount,:) = ...
            [ ...
            currentA, ...
            nextA, ...
            currentB ...
            ];


        stepA = ...
            stepA+1;



    %% --------------------------------------------------------------------
    % B 推进
    %% --------------------------------------------------------------------

    else

        nextB = ...
            ringB( ...
                mod(stepB+1,nB)+1);


        faceCount = ...
            faceCount+1;


        faces(faceCount,:) = ...
            [ ...
            currentA, ...
            nextB, ...
            currentB ...
            ];


        stepB = ...
            stepB+1;

    end

end


faces = ...
    faces(1:faceCount,:);

end



%% ========================================================================
% LOCAL FUNCTION
% 封闭端部
%% ========================================================================
function [V, ...
          faces, ...
          capNodeID, ...
          capType, ...
          nodeSectionID, ...
          nodeCircumID, ...
          nodeRingFraction, ...
          nodeIsCap] = ...
    closeBodyEnd( ...
        V, ...
        rawMesh, ...
        ringNodeIDs, ...
        side, ...
        p, ...
        nodeSectionID, ...
        nodeCircumID, ...
        nodeRingFraction, ...
        nodeIsCap)

Vraw = ...
    double(rawMesh.V);


ringNodeIDs = ...
    ringNodeIDs(:);


ringPoints = ...
    V(ringNodeIDs,:);


ringCenterYZ = ...
    mean( ...
        ringPoints(:,2:3), ...
        1);


ringRadius = ...
    mean( ...
        sqrt( ...
            (ringPoints(:,2)-ringCenterYZ(1)).^2 ...
            + ...
            (ringPoints(:,3)-ringCenterYZ(2)).^2));


xMin = ...
    min(Vraw(:,1));


xMax = ...
    max(Vraw(:,1));


bodyLength = ...
    xMax-xMin;


extremeTolerance = ...
    max( ...
        p.extremeVertexToleranceFraction ...
        * ...
        bodyLength, ...
        100 ...
        * ...
        eps(bodyLength));


switch lower(side)

    case 'aft'

        xExtreme = ...
            xMin;


        extremeMask = ...
            abs(Vraw(:,1)-xMin) ...
            <= ...
            extremeTolerance;


    case 'nose'

        xExtreme = ...
            xMax;


        extremeMask = ...
            abs(Vraw(:,1)-xMax) ...
            <= ...
            extremeTolerance;


    otherwise

        error( ...
            'RocketAeroHeat:Remesh:UnknownEnd', ...
            'side 必须为 aft 或 nose。');

end


extremeVertices = ...
    Vraw(extremeMask,:);


if isempty(extremeVertices)

    [~,index] = ...
        min( ...
            abs(Vraw(:,1)-xExtreme));


    extremeVertices = ...
        Vraw(index,:);

end


extremeCenterYZ = ...
    mean( ...
        extremeVertices(:,2:3), ...
        1);


extremeRadius = ...
    max( ...
        sqrt( ...
            (extremeVertices(:,2)-extremeCenterYZ(1)).^2 ...
            + ...
            (extremeVertices(:,3)-extremeCenterYZ(2)).^2));


if isempty(extremeRadius) ...
        || ...
        ~isfinite(extremeRadius)

    extremeRadius = ...
        0;

end


radiusRatio = ...
    extremeRadius ...
    / ...
    max(ringRadius,eps);



%% ------------------------------------------------------------------------
% 平端
%% ------------------------------------------------------------------------

if radiusRatio >= p.flatEndRadiusRatio

    capPoint = ...
        [ ...
        xExtreme, ...
        extremeCenterYZ(1), ...
        extremeCenterYZ(2) ...
        ];


    capType = ...
        'flat';



%% ------------------------------------------------------------------------
% 尖端
%% ------------------------------------------------------------------------

else

    capPoint = ...
        mean( ...
            extremeVertices, ...
            1);


    capPoint = ...
        reshape( ...
            capPoint, ...
            1, ...
            3);


    capPoint(1) = ...
        xExtreme;


    capType = ...
        'tip';

end



%% ------------------------------------------------------------------------
% 新 cap node
%% ------------------------------------------------------------------------

capNodeID = ...
    size(V,1)+1;


V(end+1,:) = ...
    reshape( ...
        capPoint, ...
        1, ...
        3);


nodeSectionID(end+1,1) = ...
    0;


nodeCircumID(end+1,1) = ...
    0;


nodeRingFraction(end+1,1) = ...
    NaN;


nodeIsCap(end+1,1) = ...
    true;


faces = ...
    makeCapFan( ...
        ringNodeIDs, ...
        capNodeID);

end



%% ========================================================================
% LOCAL FUNCTION
% Cap fan
%% ========================================================================
function faces = ...
    makeCapFan(ringNodeIDs,capNodeID)

ringNodeIDs = ...
    ringNodeIDs(:);


numberPoint = ...
    numel(ringNodeIDs);


faces = ...
    zeros(numberPoint,3);


for i = 1:numberPoint

    j = ...
        mod(i,numberPoint)+1;


    faces(i,:) = ...
        [ ...
        ringNodeIDs(i), ...
        ringNodeIDs(j), ...
        capNodeID ...
        ];

end

end



%% ========================================================================
% LOCAL FUNCTION
% 初步外法向
%% ========================================================================
function F = ...
    orientFacesOutward( ...
        V, ...
        F, ...
        faceBandID, ...
        sectionX, ...
        sectionCenterYZ, ...
        numberSections)

P1 = ...
    V(F(:,1),:);


P2 = ...
    V(F(:,2),:);


P3 = ...
    V(F(:,3),:);


normal = ...
    cross( ...
        P2-P1, ...
        P3-P1, ...
        2);


faceCenter = ...
    (P1+P2+P3)/3;


flipMask = ...
    false(size(F,1),1);



%% ------------------------------------------------------------------------
% Aft
%% ------------------------------------------------------------------------

aftMask = ...
    faceBandID == -1;


if any(aftMask)

    flipMask(aftMask) = ...
        normal(aftMask,1) ...
        > ...
        0;

end



%% ------------------------------------------------------------------------
% Nose
%% ------------------------------------------------------------------------

noseMask = ...
    faceBandID == numberSections;


if any(noseMask)

    flipMask(noseMask) = ...
        normal(noseMask,1) ...
        < ...
        0;

end



%% ------------------------------------------------------------------------
% Side
%% ------------------------------------------------------------------------

sideMask = ...
    ~(aftMask | noseMask);


if any(sideMask)

    sideIndex = ...
        find(sideMask);


    x = ...
        faceCenter(sideIndex,1);


    centerY = ...
        interp1( ...
            sectionX, ...
            sectionCenterYZ(:,1), ...
            x, ...
            'linear', ...
            'extrap');


    centerZ = ...
        interp1( ...
            sectionX, ...
            sectionCenterYZ(:,2), ...
            x, ...
            'linear', ...
            'extrap');


    outwardY = ...
        faceCenter(sideIndex,2) ...
        - ...
        centerY;


    outwardZ = ...
        faceCenter(sideIndex,3) ...
        - ...
        centerZ;


    score = ...
        normal(sideIndex,2).*outwardY ...
        + ...
        normal(sideIndex,3).*outwardZ;


    flipMask(sideIndex) = ...
        score < 0;

end


F(flipMask,[2 3]) = ...
    F(flipMask,[3 2]);

end



%% ========================================================================
% LOCAL FUNCTION
% 清理网格
%% ========================================================================
function [V,F,param] = ...
    cleanGeneratedMesh( ...
        V, ...
        F, ...
        param, ...
        referenceLength)

P1 = ...
    V(F(:,1),:);


P2 = ...
    V(F(:,2),:);


P3 = ...
    V(F(:,3),:);


area2 = ...
    sqrt( ...
        sum( ...
            cross( ...
                P2-P1, ...
                P3-P1, ...
                2).^2, ...
            2));


areaTolerance2 = ...
    2e-14 ...
    * ...
    referenceLength^2;


validFace = ...
    isfinite(area2) ...
    & ...
    area2 ...
    > ...
    areaTolerance2;


F = ...
    F(validFace,:);


param.faceBandID = ...
    param.faceBandID(validFace);



%% ------------------------------------------------------------------------
% 重复三角面
%% ------------------------------------------------------------------------

if isempty(F)

    error( ...
        'RocketAeroHeat:Remesh:AllFacesRemoved', ...
        '清理退化面后没有剩余三角形。');

end


sortedFace = ...
    sort(F,2);


[~,uniqueFaceIndex] = ...
    unique( ...
        sortedFace, ...
        'rows', ...
        'stable');


F = ...
    F(uniqueFaceIndex,:);


param.faceBandID = ...
    param.faceBandID(uniqueFaceIndex);



%% ------------------------------------------------------------------------
% 删除未使用节点
%% ------------------------------------------------------------------------

usedNode = ...
    false(size(V,1),1);


usedNode(F(:)) = ...
    true;


oldToNew = ...
    zeros(size(V,1),1);


oldToNew(usedNode) = ...
    (1:sum(usedNode))';


V = ...
    V(usedNode,:);


F = ...
    oldToNew(F);


param.nodeSectionID = ...
    param.nodeSectionID(usedNode);


param.nodeCircumID = ...
    param.nodeCircumID(usedNode);


param.nodeRingFraction = ...
    param.nodeRingFraction(usedNode);


param.nodeIsCap = ...
    param.nodeIsCap(usedNode);



%% ------------------------------------------------------------------------
% Cap ID
%% ------------------------------------------------------------------------

if param.aftCapNodeID > 0 ...
        && ...
        param.aftCapNodeID <= numel(oldToNew)

    newID = ...
        oldToNew(param.aftCapNodeID);


    param.aftCapNodeID = ...
        newID;

end


if param.noseCapNodeID > 0 ...
        && ...
        param.noseCapNodeID <= numel(oldToNew)

    newID = ...
        oldToNew(param.noseCapNodeID);


    param.noseCapNodeID = ...
        newID;

end

end



%% ========================================================================
% LOCAL FUNCTION
% 拓扑检查
%% ========================================================================
function topology = ...
    checkSurfaceTopology(F)

edges = [
    F(:,[1 2]);
    F(:,[2 3]);
    F(:,[3 1])
];


edges = ...
    sort(edges,2);


[~,~,group] = ...
    unique( ...
        edges, ...
        'rows');


count = ...
    accumarray( ...
        group, ...
        1);


topology = struct;


topology.boundaryEdges = ...
    sum(count == 1);


topology.manifoldEdges = ...
    sum(count == 2);


topology.nonManifoldEdges = ...
    sum(count > 2);


topology.isClosed = ...
    topology.boundaryEdges == 0 ...
    && ...
    topology.nonManifoldEdges == 0;

end



%% ========================================================================
% LOCAL FUNCTION
% 三角形质量
%% ========================================================================
function quality = ...
    calculateTriangleQuality(V,F)

P1 = ...
    V(F(:,1),:);


P2 = ...
    V(F(:,2),:);


P3 = ...
    V(F(:,3),:);


a = ...
    sqrt( ...
        sum( ...
            (P2-P3).^2, ...
            2));


b = ...
    sqrt( ...
        sum( ...
            (P3-P1).^2, ...
            2));


c = ...
    sqrt( ...
        sum( ...
            (P1-P2).^2, ...
            2));


crossVector = ...
    cross( ...
        P2-P1, ...
        P3-P1, ...
        2);


area = ...
    0.5 ...
    * ...
    sqrt( ...
        sum( ...
            crossVector.^2, ...
            2));


% Match the public checkMeshQuality definition used by the GUI and main.
inRadius = 2 .* area ./ max(a+b+c,eps);
aspectRatio = max([a b c],[],2) ./ max(2*sqrt(3).*inRadius,eps);


cosA = ...
    ( ...
    b.^2 ...
    + ...
    c.^2 ...
    - ...
    a.^2 ...
    ) ...
    ./ ...
    max( ...
        2.*b.*c, ...
        eps);


cosB = ...
    ( ...
    c.^2 ...
    + ...
    a.^2 ...
    - ...
    b.^2 ...
    ) ...
    ./ ...
    max( ...
        2.*c.*a, ...
        eps);


cosA = ...
    max( ...
        -1, ...
        min(1,cosA));


cosB = ...
    max( ...
        -1, ...
        min(1,cosB));


angleA = ...
    acosd(cosA);


angleB = ...
    acosd(cosB);


angleC = ...
    180 ...
    - ...
    angleA ...
    - ...
    angleB;


minimumAngle = ...
    min( ...
        [ ...
        angleA, ...
        angleB, ...
        angleC ...
        ], ...
        [], ...
        2);


quality = struct;


quality.aspectRatio = ...
    aspectRatio;


quality.minAngle = ...
    minimumAngle;


quality.area = ...
    area;

end



%% ========================================================================
% LOCAL FUNCTION
% 面元几何
%% ========================================================================
function mesh = ...
    calculateBasicFaceGeometry(mesh)

V = ...
    mesh.V;


F = ...
    mesh.F;


P1 = ...
    V(F(:,1),:);


P2 = ...
    V(F(:,2),:);


P3 = ...
    V(F(:,3),:);


rawNormal = ...
    cross( ...
        P2-P1, ...
        P3-P1, ...
        2);


normalMagnitude = ...
    sqrt( ...
        sum( ...
            rawNormal.^2, ...
            2));


mesh.area = ...
    0.5 ...
    * ...
    normalMagnitude;


mesh.normal = ...
    rawNormal ...
    ./ ...
    max( ...
        normalMagnitude, ...
        eps);


mesh.center = ...
    ( ...
    P1 ...
    + ...
    P2 ...
    + ...
    P3 ...
    ) ...
    / ...
    3;


mesh.totalArea = ...
    sum(mesh.area);


mesh.boundingBox = [
    min(V(:,1)), max(V(:,1));
    min(V(:,2)), max(V(:,2));
    min(V(:,3)), max(V(:,3))
];


mesh.size = ...
    mesh.boundingBox(:,2)' ...
    - ...
    mesh.boundingBox(:,1)';


characteristicLength = ...
    max(mesh.size);


if ~(isfinite(characteristicLength) ...
        && ...
        characteristicLength > 0)

    characteristicLength = ...
        1;

end


mesh.badFace = ...
    mesh.area ...
    <= ...
    1e-14 ...
    * ...
    characteristicLength^2;


mesh.numberBadFaces = ...
    sum(mesh.badFace);

end



%% ========================================================================
% LOCAL FUNCTION
% Mesh array validation
%% ========================================================================
function validateXYZMeshArrays(V,F,label)

if nargin < 3

    label = ...
        'mesh';

end


if ~isnumeric(V) ...
        || ...
        ndims(V) ~= 2 ...
        || ...
        size(V,2) ~= 3

    error( ...
        'RocketAeroHeat:Remesh:VertexDimension', ...
        ['%s.V 必须为 N×3。' newline ...
         '当前尺寸 = %d×%d。'], ...
        label, ...
        size(V,1), ...
        size(V,2));

end


if ~isnumeric(F) ...
        || ...
        ndims(F) ~= 2 ...
        || ...
        size(F,2) ~= 3

    error( ...
        'RocketAeroHeat:Remesh:FaceDimension', ...
        ['%s.F 必须为 M×3。' newline ...
         '当前尺寸 = %d×%d。'], ...
        label, ...
        size(F,1), ...
        size(F,2));

end


if isempty(V)

    error( ...
        'RocketAeroHeat:Remesh:EmptyVertex', ...
        '%s.V 为空。', ...
        label);

end


if isempty(F)

    error( ...
        'RocketAeroHeat:Remesh:EmptyFace', ...
        '%s.F 为空。', ...
        label);

end


if any(F(:) < 1) ...
        || ...
        any(F(:) > size(V,1)) ...
        || ...
        any(F(:) ~= round(F(:)))

    error( ...
        'RocketAeroHeat:Remesh:FaceIndex', ...
        '%s.F 含非法节点编号。', ...
        label);

end

end
