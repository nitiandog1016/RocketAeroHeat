function mesh = buildSurfaceCoordinates(mesh, cfg)
%% ========================================================================
% buildSurfaceCoordinates.m
%
% Rocket AeroHeat
%
% 三维箭体表面参数坐标构建模块
%
% ========================================================================
%
% 功能：
%
% 给最终计算表面上的每一个节点、每一个三角面元建立：
%
%       x
%       r
%       theta
%       s
%
%
% ========================================================================
%
% 软件内部坐标：
%
%       +X = 鼻尖 / 前向
%       +Y = 横向
%       +Z = 向上
%
%
% ========================================================================
%
% 定义：
%
% x
%
%       节点或面元中心的轴向坐标
%
%
% r
%
%       相对于局部中心线：
%
%           r = sqrt[(y-yc)^2 + (z-zc)^2]
%
%
% theta
%
%       Y-Z 平面周向角：
%
%           theta = atan2(z-zc, y-yc)
%
%       并统一到：
%
%           0 <= theta < 2*pi
%
%
% sFromAft
%
%       从尾端向鼻端沿几何子午方向累计的表面距离
%
%
% sFromNose
%
%       从鼻端向尾端沿几何子午方向累计的表面距离
%
%
% nodeS / faceS
%
%       当前默认：
%
%           = sFromNose
%
%
% ========================================================================
%
% 非常重要：
%
% 这里的 s 是：
%
%       几何参数子午线弧长
%
% 而不是：
%
%       有攻角情况下的真实表面流线长度
%
%
% 对：
%
%       轴对称体
%       小攻角长细体
%
% 它非常有用。
%
% 后续 classify / boundary-layer 模块还会建立真正的：
%
%       streamline distance
%
%
% ========================================================================
%
% 推荐调用：
%
%       mesh = buildSurfaceCoordinates(mesh,cfg);
%
%
% 输出：
%
% mesh.surfaceCoord.nodeX
% mesh.surfaceCoord.nodeR
% mesh.surfaceCoord.nodeTheta
% mesh.surfaceCoord.nodeS
% mesh.surfaceCoord.nodeSFromNose
% mesh.surfaceCoord.nodeSFromAft
%
% mesh.surfaceCoord.faceX
% mesh.surfaceCoord.faceR
% mesh.surfaceCoord.faceTheta
% mesh.surfaceCoord.faceS
% mesh.surfaceCoord.faceSFromNose
% mesh.surfaceCoord.faceSFromAft
%
% mesh.surfaceCoord.nodeSectionID
% mesh.surfaceCoord.nodeCircumID
% mesh.surfaceCoord.nodeRingFraction
%
% mesh.surfaceCoord.sectionX
% mesh.surfaceCoord.sectionCenterYZ
%
% mesh.surfaceCoord.commonRingFraction
% mesh.surfaceCoord.sectionSFromNoseGrid
% mesh.surfaceCoord.sectionSFromAftGrid
%
%% ========================================================================



%% ========================================================================
% 1. 输入
% ========================================================================

if nargin < 1

    error( ...
        'buildSurfaceCoordinates: 必须输入 mesh。');

end


if nargin < 2

    cfg = struct;

end


if ~isfield(mesh,'V') ...
        || ...
        ~isfield(mesh,'F')

    error( ...
        'buildSurfaceCoordinates: mesh 必须包含 V 和 F。');

end


V = ...
    double(mesh.V);


F = ...
    double(mesh.F);


numberNode = ...
    size(V,1);


numberFace = ...
    size(F,1);



%% ========================================================================
% 2. 检查 / 获取半结构化拓扑信息
% ========================================================================

hasNativeTopology = ...
    isfield(mesh,'param') ...
    && ...
    isfield(mesh.param,'sectionNodeIDs') ...
    && ...
    isfield(mesh.param,'sectionX') ...
    && ...
    isfield(mesh.param,'nodeSectionID') ...
    && ...
    isfield(mesh.param,'nodeRingFraction');


if hasNativeTopology

    param = ...
        mesh.param;


else

    warning( ...
        ['buildSurfaceCoordinates: mesh.param 不完整。' newline ...
         '将根据 X 坐标尝试反推轴向截面。' newline ...
         '建议使用新版 autoRemeshSurface.m 重新生成网格。']);


    param = ...
        inferSectionTopologyFromMesh( ...
            mesh, ...
            cfg);

end



%% ========================================================================
% 3. 读取截面信息
% ========================================================================

sectionX = ...
    param.sectionX(:);


sectionNodeIDs = ...
    param.sectionNodeIDs;


numberSection = ...
    numel(sectionX);


if numberSection < 2

    error( ...
        'buildSurfaceCoordinates: 有效截面少于 2。');

end



%% ========================================================================
% 4. 截面中心线
% ========================================================================

if isfield(param,'sectionCenterYZ') ...
        && ...
        size(param.sectionCenterYZ,1) == numberSection

    sectionCenterYZ = ...
        param.sectionCenterYZ;


else

    sectionCenterYZ = ...
        zeros(numberSection,2);


    for i = 1:numberSection

        ids = ...
            sectionNodeIDs{i};


        sectionCenterYZ(i,:) = ...
            mean( ...
                V(ids,2:3), ...
                1);

    end

end



%% ========================================================================
% 5. 节点 x
% ========================================================================

nodeX = ...
    V(:,1);



%% ========================================================================
% 6. 节点局部中心线坐标
% ========================================================================

centerYNode = ...
    interp1( ...
        sectionX, ...
        sectionCenterYZ(:,1), ...
        nodeX, ...
        'linear', ...
        'extrap');


centerZNode = ...
    interp1( ...
        sectionX, ...
        sectionCenterYZ(:,2), ...
        nodeX, ...
        'linear', ...
        'extrap');



%% ========================================================================
% 7. 节点半径
% ========================================================================

nodeYRelative = ...
    V(:,2) ...
    - ...
    centerYNode;


nodeZRelative = ...
    V(:,3) ...
    - ...
    centerZNode;


nodeR = ...
    sqrt( ...
        nodeYRelative.^2 ...
        + ...
        nodeZRelative.^2);



%% ========================================================================
% 8. 节点 theta
%
% theta = 0：
%
%       +Y
%
% theta = pi/2：
%
%       +Z
%% ========================================================================

nodeTheta = ...
    atan2( ...
        nodeZRelative, ...
        nodeYRelative);


nodeTheta( ...
    nodeTheta < 0) = ...
    nodeTheta( ...
        nodeTheta < 0) ...
    + ...
    2*pi;



%% ========================================================================
% 9. 节点 topology
% ========================================================================

nodeSectionID = ...
    param.nodeSectionID(:);


if isfield(param,'nodeCircumID')

    nodeCircumID = ...
        param.nodeCircumID(:);

else

    nodeCircumID = ...
        zeros(numberNode,1);

end


nodeRingFraction = ...
    param.nodeRingFraction(:);



%% ========================================================================
% 10. 每个截面实际节点数
% ========================================================================

sectionPointCount = ...
    zeros(numberSection,1);


for i = 1:numberSection

    sectionPointCount(i) = ...
        numel(sectionNodeIDs{i});

end



%% ========================================================================
% 11. 建立公共周向参数网格
%
% 不同截面的 Ntheta 可以不同。
%
% 为了计算沿相同参数线发展的 meridional s：
%
% 先将所有截面映射到统一 u 网格：
%
%       u = theta-like ring fraction
%
%       0 <= u < 1
%% ========================================================================

numberCommonPoint = ...
    max(sectionPointCount);


numberCommonPoint = ...
    max( ...
        numberCommonPoint, ...
        8);


commonRingFraction = ...
    ( ...
    0:numberCommonPoint-1 ...
    )' ...
    / ...
    numberCommonPoint;



%% ========================================================================
% 12. 将每个实际截面插值到公共 u 网格
%
% sectionPointGrid：
%
%       numberSection × numberCommonPoint × 3
%% ========================================================================

sectionPointGrid = ...
    zeros( ...
        numberSection, ...
        numberCommonPoint, ...
        3);


for i = 1:numberSection

    ids = ...
        sectionNodeIDs{i};


    if isempty(ids)

        error( ...
            '第 %d 个截面没有节点。', ...
            i);

    end


    P = ...
        V(ids,:);


    if numel(ids) >= 2 ...
            && ...
            all(isfinite(nodeRingFraction(ids)))

        u = ...
            nodeRingFraction(ids);


    else

        u = ...
            (0:numel(ids)-1)' ...
            / ...
            numel(ids);

    end


    %% --------------------------------------------------------------------
    % 确保按 u 排序
    %% --------------------------------------------------------------------

    [u,order] = ...
        sort(u);


    P = ...
        P(order,:);


    for coordinate = 1:3

        sectionPointGrid(i,:,coordinate) = ...
            periodicInterpolate( ...
                u, ...
                P(:,coordinate), ...
                commonRingFraction);

    end

end



%% ========================================================================
% 13. 检查尾端 / 鼻端 cap node
% ========================================================================

aftCapPoint = [];

noseCapPoint = [];


if isfield(param,'aftCapNodeID') ...
        && ...
        param.aftCapNodeID > 0 ...
        && ...
        param.aftCapNodeID <= numberNode

    aftCapPoint = ...
        V(param.aftCapNodeID,:);

end


if isfield(param,'noseCapNodeID') ...
        && ...
        param.noseCapNodeID > 0 ...
        && ...
        param.noseCapNodeID <= numberNode

    noseCapPoint = ...
        V(param.noseCapNodeID,:);

end



%% ========================================================================
% 14. 环与环之间的 meridional distance
%
% 对每一个公共 u：
%
%       section 1
%           ↓
%       section 2
%           ↓
%       ...
%           ↓
%       section N
%% ========================================================================

sectionDistance = ...
    zeros( ...
        numberSection-1, ...
        numberCommonPoint);


for i = 1:numberSection-1

    P1 = ...
        squeeze( ...
            sectionPointGrid(i,:,:));


    P2 = ...
        squeeze( ...
            sectionPointGrid(i+1,:,:));


    sectionDistance(i,:) = ...
        sqrt( ...
            sum( ...
                (P2-P1).^2, ...
                2))';

end



%% ========================================================================
% 15. 尾端延伸距离
%% ========================================================================

aftExtension = ...
    zeros(1,numberCommonPoint);


if ~isempty(aftCapPoint)

    firstRing = ...
        squeeze( ...
            sectionPointGrid(1,:,:));


    aftExtension = ...
        sqrt( ...
            sum( ...
                ( ...
                firstRing ...
                - ...
                aftCapPoint ...
                ).^2, ...
                2))';

end



%% ========================================================================
% 16. 鼻端延伸距离
%% ========================================================================

noseExtension = ...
    zeros(1,numberCommonPoint);


if ~isempty(noseCapPoint)

    lastRing = ...
        squeeze( ...
            sectionPointGrid(end,:,:));


    noseExtension = ...
        sqrt( ...
            sum( ...
                ( ...
                lastRing ...
                - ...
                noseCapPoint ...
                ).^2, ...
                2))';

end



%% ========================================================================
% 17. sFromAft 网格
%% ========================================================================

sectionSFromAftGrid = ...
    zeros( ...
        numberSection, ...
        numberCommonPoint);


sectionSFromAftGrid(1,:) = ...
    aftExtension;


for i = 2:numberSection

    sectionSFromAftGrid(i,:) = ...
        sectionSFromAftGrid(i-1,:) ...
        + ...
        sectionDistance(i-1,:);

end



%% ========================================================================
% 18. 每条参数子午线总长度
%% ========================================================================

totalMeridionalLength = ...
    sectionSFromAftGrid(end,:) ...
    + ...
    noseExtension;



%% ========================================================================
% 19. 从鼻端量起
%% ========================================================================

sectionSFromNoseGrid = ...
    zeros( ...
        numberSection, ...
        numberCommonPoint);


for i = 1:numberSection

    sectionSFromNoseGrid(i,:) = ...
        totalMeridionalLength ...
        - ...
        sectionSFromAftGrid(i,:);

end



%% ========================================================================
% 20. 将公共 s 网格重新插值回原始节点
%% ========================================================================

nodeSFromAft = ...
    nan(numberNode,1);


nodeSFromNose = ...
    nan(numberNode,1);


for i = 1:numberSection

    ids = ...
        sectionNodeIDs{i};


    if isempty(ids)

        continue;

    end


    uNode = ...
        nodeRingFraction(ids);


    invalidU = ...
        ~isfinite(uNode);


    if any(invalidU)

        numberLocal = ...
            numel(ids);


        uNode(invalidU) = ...
            ( ...
            find(invalidU)-1 ...
            ) ...
            / ...
            numberLocal;

    end


    nodeSFromAft(ids) = ...
        periodicInterpolate( ...
            commonRingFraction, ...
            sectionSFromAftGrid(i,:)', ...
            uNode);


    nodeSFromNose(ids) = ...
        periodicInterpolate( ...
            commonRingFraction, ...
            sectionSFromNoseGrid(i,:)', ...
            uNode);

end



%% ========================================================================
% 21. Cap 节点 s
%% ========================================================================

if ~isempty(aftCapPoint) ...
        && ...
        isfield(param,'aftCapNodeID') ...
        && ...
        param.aftCapNodeID > 0


    aftID = ...
        param.aftCapNodeID;


    nodeSFromAft(aftID) = ...
        0;


    nodeSFromNose(aftID) = ...
        mean(totalMeridionalLength);

end


if ~isempty(noseCapPoint) ...
        && ...
        isfield(param,'noseCapNodeID') ...
        && ...
        param.noseCapNodeID > 0


    noseID = ...
        param.noseCapNodeID;


    nodeSFromNose(noseID) = ...
        0;


    nodeSFromAft(noseID) = ...
        mean(totalMeridionalLength);

end



%% ========================================================================
% 22. 未定义节点备用处理
%% ========================================================================

undefinedAft = ...
    ~isfinite(nodeSFromAft);


undefinedNose = ...
    ~isfinite(nodeSFromNose);


if any(undefinedAft) ...
        || ...
        any(undefinedNose)


    xMinimum = ...
        min(nodeX);


    xMaximum = ...
        max(nodeX);


    normalizedX = ...
        ( ...
        nodeX-xMinimum ...
        ) ...
        / ...
        max(xMaximum-xMinimum,eps);


    meanLength = ...
        mean(totalMeridionalLength);


    nodeSFromAft(undefinedAft) = ...
        normalizedX(undefinedAft) ...
        * ...
        meanLength;


    nodeSFromNose(undefinedNose) = ...
        ( ...
        1-normalizedX(undefinedNose) ...
        ) ...
        * ...
        meanLength;

end



%% ========================================================================
% 23. 当前默认 s
%
% 因为：
%
%       +X = nose
%
% 后续零攻角高超声速边界层通常从鼻部开始发展。
%% ========================================================================

nodeS = ...
    nodeSFromNose;



%% ========================================================================
% 24. 面元中心
%% ========================================================================

if isfield(mesh,'center') ...
        && ...
        size(mesh.center,1) == numberFace

    faceCenter = ...
        mesh.center;


else

    faceCenter = ...
        ( ...
        V(F(:,1),:) ...
        + ...
        V(F(:,2),:) ...
        + ...
        V(F(:,3),:) ...
        ) ...
        / ...
        3;

end



%% ========================================================================
% 25. 面元 X
%% ========================================================================

faceX = ...
    faceCenter(:,1);



%% ========================================================================
% 26. 面元局部中心线
%% ========================================================================

centerYFace = ...
    interp1( ...
        sectionX, ...
        sectionCenterYZ(:,1), ...
        faceX, ...
        'linear', ...
        'extrap');


centerZFace = ...
    interp1( ...
        sectionX, ...
        sectionCenterYZ(:,2), ...
        faceX, ...
        'linear', ...
        'extrap');



%% ========================================================================
% 27. 面元 r
%% ========================================================================

faceYRelative = ...
    faceCenter(:,2) ...
    - ...
    centerYFace;


faceZRelative = ...
    faceCenter(:,3) ...
    - ...
    centerZFace;


faceR = ...
    sqrt( ...
        faceYRelative.^2 ...
        + ...
        faceZRelative.^2);



%% ========================================================================
% 28. 面元 theta
%% ========================================================================

faceTheta = ...
    atan2( ...
        faceZRelative, ...
        faceYRelative);


faceTheta( ...
    faceTheta < 0) = ...
    faceTheta( ...
        faceTheta < 0) ...
    + ...
    2*pi;



%% ========================================================================
% 29. 面元 s
%
% 使用三个顶点 s 的平均值。
%% ========================================================================

faceSFromAft = ...
    ( ...
    nodeSFromAft(F(:,1)) ...
    + ...
    nodeSFromAft(F(:,2)) ...
    + ...
    nodeSFromAft(F(:,3)) ...
    ) ...
    / ...
    3;


faceSFromNose = ...
    ( ...
    nodeSFromNose(F(:,1)) ...
    + ...
    nodeSFromNose(F(:,2)) ...
    + ...
    nodeSFromNose(F(:,3)) ...
    ) ...
    / ...
    3;


faceS = ...
    faceSFromNose;



%% ========================================================================
% 30. 面元 Section / Band ID
%% ========================================================================

if isfield(param,'faceBandID') ...
        && ...
        numel(param.faceBandID) == numberFace

    faceBandID = ...
        param.faceBandID(:);


else

    faceBandID = ...
        estimateFaceBandID( ...
            F, ...
            nodeSectionID, ...
            numberSection);

end



%% ========================================================================
% 31. 面元周向参数
%
% 使用上面按面心与局部中心线计算的 atan2 方位角。
% 顶点角的单位圆平均会给轴线上没有方位角的鼻尖赋予额外权重，
% 并且不等于不同半径顶点所组成三角形的面心方位角。
%% ========================================================================




%% ========================================================================
% 32. 保存结果
%% ========================================================================

surfaceCoord = struct;


surfaceCoord.version = ...
    1;


surfaceCoord.description = ...
    ['Body-fitted axial/circumferential/meridional ', ...
     'surface coordinate system'];


surfaceCoord.hasNativeSectionTopology = ...
    hasNativeTopology;



%% ------------------------------------------------------------------------
% 坐标约定
%% ------------------------------------------------------------------------

surfaceCoord.bodyAxis = ...
    [1 0 0];


surfaceCoord.noseDirection = ...
    [1 0 0];


surfaceCoord.thetaDefinition = ...
    'theta = atan2(z-zc, y-yc), range [0,2*pi)';


surfaceCoord.sDefinition = ...
    ['Geometrical meridional distance based on constant ', ...
     'ring-fraction parameter lines. Not actual streamline distance.'];



%% ------------------------------------------------------------------------
% Node
%% ------------------------------------------------------------------------

surfaceCoord.nodeX = ...
    nodeX;


surfaceCoord.nodeR = ...
    nodeR;


surfaceCoord.nodeTheta = ...
    nodeTheta;


surfaceCoord.nodeS = ...
    nodeS;


surfaceCoord.nodeSFromNose = ...
    nodeSFromNose;


surfaceCoord.nodeSFromAft = ...
    nodeSFromAft;


surfaceCoord.nodeSectionID = ...
    nodeSectionID;


surfaceCoord.nodeCircumID = ...
    nodeCircumID;


surfaceCoord.nodeRingFraction = ...
    nodeRingFraction;


if isfield(param,'nodeIsCap')

    surfaceCoord.nodeIsCap = ...
        param.nodeIsCap;

else

    surfaceCoord.nodeIsCap = ...
        nodeSectionID == 0;

end



%% ------------------------------------------------------------------------
% Face
%% ------------------------------------------------------------------------

surfaceCoord.faceX = ...
    faceX;


surfaceCoord.faceR = ...
    faceR;


surfaceCoord.faceTheta = ...
    faceTheta;


surfaceCoord.faceS = ...
    faceS;


surfaceCoord.faceSFromNose = ...
    faceSFromNose;


surfaceCoord.faceSFromAft = ...
    faceSFromAft;


surfaceCoord.faceBandID = ...
    faceBandID;



%% ------------------------------------------------------------------------
% Section
%% ------------------------------------------------------------------------

surfaceCoord.sectionX = ...
    sectionX;


surfaceCoord.sectionCenterYZ = ...
    sectionCenterYZ;


surfaceCoord.sectionPointCount = ...
    sectionPointCount;


surfaceCoord.sectionNodeIDs = ...
    sectionNodeIDs;



%% ------------------------------------------------------------------------
% 公共参数网格
%% ------------------------------------------------------------------------

surfaceCoord.commonRingFraction = ...
    commonRingFraction;


surfaceCoord.sectionPointGrid = ...
    sectionPointGrid;


surfaceCoord.sectionSFromAftGrid = ...
    sectionSFromAftGrid;


surfaceCoord.sectionSFromNoseGrid = ...
    sectionSFromNoseGrid;


surfaceCoord.totalMeridionalLengthByFraction = ...
    totalMeridionalLength(:);



%% ------------------------------------------------------------------------
% 标量统计
%% ------------------------------------------------------------------------

surfaceCoord.meanMeridionalLength = ...
    mean(totalMeridionalLength);


surfaceCoord.minimumMeridionalLength = ...
    min(totalMeridionalLength);


surfaceCoord.maximumMeridionalLength = ...
    max(totalMeridionalLength);


surfaceCoord.numberSections = ...
    numberSection;


surfaceCoord.numberCommonCircumPoints = ...
    numberCommonPoint;



%% ========================================================================
% 33. 写回 mesh
%% ========================================================================

mesh.surfaceCoord = ...
    surfaceCoord;



%% ========================================================================
% 34. 输出报告
%% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf(' Surface Coordinate Construction\n');
fprintf('============================================================\n');


fprintf('Native section topology  : %s\n', ...
    logicalText(hasNativeTopology));


fprintf('Axial sections           : %d\n', ...
    numberSection);


fprintf('Common circum points     : %d\n', ...
    numberCommonPoint);


fprintf('\n');


fprintf('Mean meridional length   : %.6f m\n', ...
    surfaceCoord.meanMeridionalLength);


fprintf('Minimum meridional length: %.6f m\n', ...
    surfaceCoord.minimumMeridionalLength);


fprintf('Maximum meridional length: %.6f m\n', ...
    surfaceCoord.maximumMeridionalLength);


fprintf('\n');


fprintf('Default s coordinate     : From nose toward aft\n');


fprintf('theta = 0                : +Y direction\n');

fprintf('theta = pi/2             : +Z direction\n');


fprintf('\n');


fprintf('NOTE:\n');

fprintf('s is a geometric meridional coordinate.\n');

fprintf('It is NOT yet the true surface streamline distance at AoA/beta.\n');


fprintf('============================================================\n');

end



%% ========================================================================
% LOCAL FUNCTION
% 周期插值
%% ========================================================================
function valueQuery = ...
    periodicInterpolate(u,value,uQuery)

u = ...
    u(:);


value = ...
    value(:);


uQuery = ...
    uQuery(:);


if numel(u) ~= numel(value)

    error( ...
        'periodicInterpolate: u 和 value 长度不一致。');

end


%% ------------------------------------------------------------------------
% 归一化到 [0,1)
%% ------------------------------------------------------------------------

u = ...
    mod(u,1);


uQuery = ...
    mod(uQuery,1);


%% ------------------------------------------------------------------------
% 排序
%% ------------------------------------------------------------------------

[u,order] = ...
    sort(u);


value = ...
    value(order);



%% ------------------------------------------------------------------------
% 删除重复 u
%% ------------------------------------------------------------------------

[u,uniqueIndex] = ...
    unique( ...
        u, ...
        'stable');


value = ...
    value(uniqueIndex);


if numel(u) == 1

    valueQuery = ...
        value(1) ...
        * ...
        ones(size(uQuery));

    return;

end



%% ------------------------------------------------------------------------
% 周期闭合
%% ------------------------------------------------------------------------

uExtended = [
    u;
    u(1)+1
];


valueExtended = [
    value;
    value(1)
];


%% ------------------------------------------------------------------------
% 如果第一个 u 不是 0：
%
% 在负方向加上最后一个周期点。
%% ------------------------------------------------------------------------

if u(1) > 1e-12

    uExtended = [
        u(end)-1;
        uExtended
    ];


    valueExtended = [
        value(end);
        valueExtended
    ];

end



%% ------------------------------------------------------------------------
% 插值
%% ------------------------------------------------------------------------

valueQuery = ...
    interp1( ...
        uExtended, ...
        valueExtended, ...
        uQuery, ...
        'linear');

end



%% ========================================================================
% LOCAL FUNCTION
% 老网格没有 param 时根据 X 反推截面
%% ========================================================================
function param = ...
    inferSectionTopologyFromMesh(mesh,cfg)

V = ...
    mesh.V;


numberNode = ...
    size(V,1);


modelRange = ...
    max(V,[],1) ...
    - ...
    min(V,[],1);


referenceLength = ...
    max(modelRange);


tolerance = ...
    getCfgValue( ...
        cfg, ...
        {'remesh','sectionMergeToleranceFraction'}, ...
        1e-7) ...
    * ...
    referenceLength;


tolerance = ...
    max( ...
        tolerance, ...
        100*eps(referenceLength));



%% ========================================================================
% X 排序
%% ========================================================================

[xSorted,order] = ...
    sort(V(:,1));


groups = ...
    zeros(numberNode,1);


groupNumber = ...
    1;


groups(order(1)) = ...
    groupNumber;


for k = 2:numberNode

    if abs( ...
            xSorted(k) ...
            - ...
            xSorted(k-1)) ...
            > ...
            tolerance

        groupNumber = ...
            groupNumber+1;

    end


    groups(order(k)) = ...
        groupNumber;

end



%% ========================================================================
% 只保留至少 3 个节点的环
%% ========================================================================

% A flat-cap center shares X with its perimeter, but is not a ring node.
faceX = reshape(V(mesh.F(:),1),size(mesh.F));
sideFaces = max(faceX,[],2)-min(faceX,[],2) > tolerance;
sideNodes = unique(mesh.F(sideFaces,:));
atEnd = abs(V(:,1)-min(V(:,1))) <= tolerance | ...
    abs(V(:,1)-max(V(:,1))) <= tolerance;
ringCandidate = ~atEnd | ismember((1:numberNode)',sideNodes);
groupCount = accumarray(groups(ringCandidate),1,[groupNumber 1]);


ringGroups = ...
    find(groupCount >= 3);


numberSection = ...
    numel(ringGroups);


sectionX = ...
    zeros(numberSection,1);


sectionCenterYZ = ...
    zeros(numberSection,2);


sectionNodeIDs = ...
    cell(numberSection,1);


nodeSectionID = ...
    zeros(numberNode,1);


nodeCircumID = ...
    zeros(numberNode,1);


nodeRingFraction = ...
    nan(numberNode,1);


for section = 1:numberSection

    groupID = ...
        ringGroups(section);


    ids = ...
        find(groups == groupID & ringCandidate);


    centerYZ = ...
        mean(V(ids,2:3),1);


    theta = ...
        atan2( ...
            V(ids,3)-centerYZ(2), ...
            V(ids,2)-centerYZ(1));


    theta(theta < 0) = ...
        theta(theta < 0) ...
        + ...
        2*pi;


    [~,circumOrder] = ...
        sort(theta);


    ids = ...
        ids(circumOrder);


    numberLocal = ...
        numel(ids);


    sectionX(section) = ...
        mean(V(ids,1));


    sectionCenterYZ(section,:) = ...
        centerYZ;


    sectionNodeIDs{section} = ...
        ids;


    nodeSectionID(ids) = ...
        section;


    nodeCircumID(ids) = ...
        (1:numberLocal)';


    nodeRingFraction(ids) = ...
        (0:numberLocal-1)' ...
        / ...
        numberLocal;

end



%% ========================================================================
% 截面按 X 排序
%% ========================================================================

[sectionX,sectionOrder] = ...
    sort(sectionX);


oldSectionNodeIDs = ...
    sectionNodeIDs;


oldCenterYZ = ...
    sectionCenterYZ;


sectionNodeIDs = ...
    cell(numberSection,1);


sectionCenterYZ = ...
    zeros(numberSection,2);


newNodeSectionID = ...
    zeros(numberNode,1);


for newID = 1:numberSection

    oldID = ...
        sectionOrder(newID);


    ids = ...
        oldSectionNodeIDs{oldID};


    sectionNodeIDs{newID} = ...
        ids;


    sectionCenterYZ(newID,:) = ...
        oldCenterYZ(oldID,:);


    newNodeSectionID(ids) = ...
        newID;

end


nodeSectionID = ...
    newNodeSectionID;



%% ========================================================================
% Cap nodes
%% ========================================================================

capNodes = ...
    find(nodeSectionID == 0);


aftCapNodeID = ...
    0;


noseCapNodeID = ...
    0;


if ~isempty(capNodes)
    aftCandidates = capNodes(abs(V(capNodes,1)-min(V(:,1))) <= tolerance);
    noseCandidates = capNodes(abs(V(capNodes,1)-max(V(:,1))) <= tolerance);
    if numel(aftCandidates) > 1 || numel(noseCandidates) > 1
        error('buildSurfaceCoordinates:AmbiguousCaps', ...
            'Cannot infer a unique cap center; provide mesh.param for this topology.');
    end
    if ~isempty(aftCandidates)
        aftCapNodeID = aftCandidates(1);
    end
    if ~isempty(noseCandidates)
        noseCapNodeID = noseCandidates(1);
    end

end



%% ========================================================================
% Face band
%% ========================================================================

faceBandID = ...
    estimateFaceBandID( ...
        mesh.F, ...
        nodeSectionID, ...
        numberSection);



%% ========================================================================
% 输出
%% ========================================================================

param = struct;


param.sectionX = ...
    sectionX;


param.sectionCenterYZ = ...
    sectionCenterYZ;


param.sectionNodeIDs = ...
    sectionNodeIDs;


param.nodeSectionID = ...
    nodeSectionID;


param.nodeCircumID = ...
    nodeCircumID;


param.nodeRingFraction = ...
    nodeRingFraction;


param.nodeIsCap = ...
    nodeSectionID == 0;


param.faceBandID = ...
    faceBandID;


param.aftCapNodeID = ...
    aftCapNodeID;


param.noseCapNodeID = ...
    noseCapNodeID;

end



%% ========================================================================
% LOCAL FUNCTION
% 根据节点 section ID 推测 face band
%% ========================================================================
function faceBandID = ...
    estimateFaceBandID(F,nodeSectionID,numberSection)

numberFace = ...
    size(F,1);


faceBandID = ...
    zeros(numberFace,1);


for i = 1:numberFace

    section = ...
        nodeSectionID(F(i,:));


    nonzeroSection = ...
        section(section > 0);


    if isempty(nonzeroSection)

        faceBandID(i) = ...
            0;

        continue;

    end


    minimumSection = ...
        min(nonzeroSection);


    maximumSection = ...
        max(nonzeroSection);


    if any(section == 0)

        if maximumSection <= 1

            faceBandID(i) = ...
                -1;

        elseif minimumSection >= numberSection

            faceBandID(i) = ...
                numberSection;

        else

            faceBandID(i) = ...
                round(mean(nonzeroSection));

        end


    else

        faceBandID(i) = ...
            minimumSection;

    end

end

end



%% ========================================================================
% LOCAL FUNCTION
% cfg
%% ========================================================================
function value = ...
    getCfgValue(cfg,path,defaultValue)

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
% logical → 字符
%% ========================================================================
function str = logicalText(value)

if value

    str = ...
        'YES';

else

    str = ...
        'NO';

end

end
