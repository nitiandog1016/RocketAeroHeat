function mesh = fixNormalDirection(mesh, cfg)
%% ========================================================================
% fixNormalDirection.m
%
% Rocket AeroHeat
% 三角表面网格法向一致性与外法向修正程序
%
% ========================================================================
%
% 功能：
%
%   1. 检查三角面之间的拓扑连接关系
%
%   2. 通过公共边判断相邻三角面的方向是否一致
%
%   3. 自动翻转方向错误的三角面
%
%   4. 对每一个独立连通区域判断整体法向是否朝外
%
%   5. 对封闭表面利用有向体积判断：
%
%           signed volume > 0
%
%       表示整体方向正确。
%
%   6. 对开放表面无法使用体积时，
%      使用“面法向相对于区域中心”的方法作为备用判断。
%
%   7. 最后重新计算：
%
%           mesh.center
%           mesh.normal
%           mesh.area
%           mesh.totalArea
%
%
% ------------------------------------------------------------------------
% 输入：
%
% mesh.V
%       N × 3 节点坐标
%
% mesh.F
%       M × 3 三角面连接关系
%
% cfg
%       总配置结构体
%
%
% ------------------------------------------------------------------------
% 输出：
%
% mesh.F
%       修正方向后的三角面
%
% mesh.normalInfo
%       法向修正统计信息
%
%
% ------------------------------------------------------------------------
% 对后续程序的重要性：
%
% Modified Newtonian：
%
%       mu = -n · Vinf
%
% 如果 n 错了：
%
%       迎风面
%
% 会被错误判断成：
%
%       背风面
%
% 因此该模块应放在正式气动力计算之前。
%
%% ========================================================================


%% ========================================================================
% 0. 配置参数
% ========================================================================

verbose = ...
    getCfgValue( ...
        cfg, ...
        {'normal','verbose'}, ...
        true);


failOnNonManifold = ...
    getCfgValue( ...
        cfg, ...
        {'normal','failOnNonManifold'}, ...
        false);


%% ========================================================================
% 1. 读取网格
% ========================================================================

V = double(mesh.V);

F = double(mesh.F);


numberVertices = ...
    size(V,1);


numberFaces = ...
    size(F,1);


if numberFaces == 0

    error('fixNormalDirection: mesh.F 为空。');

end


if any(F(:) < 1) || ...
        any(F(:) > numberVertices)

    error('fixNormalDirection: mesh.F 中存在非法节点编号。');

end


%% ========================================================================
% 2. 检查一个三角面内部是否存在重复节点
% ========================================================================

invalidFace = ...
    F(:,1) == F(:,2) ...
    | ...
    F(:,2) == F(:,3) ...
    | ...
    F(:,3) == F(:,1);


if any(invalidFace)

    error( ...
        ['fixNormalDirection: 网格中存在 ', ...
         '%d 个节点编号重复的退化三角面。'], ...
         sum(invalidFace));

end


%% ========================================================================
% 3. 建立所有有向边
% ========================================================================
%
% 对三角形：
%
%       [a b c]
%
% 三条有向边为：
%
%       a → b
%       b → c
%       c → a
%
%
% 对一个方向一致的封闭三角网格：
%
% 相邻两个三角形在公共边上的方向应该相反。
%

directedEdges = [
    F(:,[1 2]);
    F(:,[2 3]);
    F(:,[3 1])
];


faceID = ...
    repmat( ...
        (1:numberFaces)', ...
        3, ...
        1);


%% ========================================================================
% 4. 将公共边转换为无向边
% ========================================================================

sortedEdges = ...
    sort( ...
        directedEdges, ...
        2);


%% ========================================================================
% 5. 保存原始边相对于排序边的方向
% ========================================================================
%
% 如果：
%
%       directed edge = 小节点 → 大节点
%
% directionSign = +1
%
%
% 如果：
%
%       directed edge = 大节点 → 小节点
%
% directionSign = -1
%

directionSign = ...
    ones( ...
        size(directedEdges,1), ...
        1);


directionSign( ...
    directedEdges(:,1) > ...
    directedEdges(:,2)) = -1;


%% ========================================================================
% 6. 对公共边进行分组
% ========================================================================

[uniqueEdges,~,edgeGroup] = ...
    unique( ...
        sortedEdges, ...
        'rows');


numberEdges = ...
    size(uniqueEdges,1);


edgeCount = ...
    accumarray( ...
        edgeGroup, ...
        1);


%% ========================================================================
% 7. 每条无向边对应哪些 edge occurrence
% ========================================================================

edgeOccurrence = ...
    accumarray( ...
        edgeGroup, ...
        (1:numel(edgeGroup))', ...
        [numberEdges 1], ...
        @(x){x});


%% ========================================================================
% 8. 网格拓扑统计
% ========================================================================

boundaryEdgeFlag = ...
    edgeCount == 1;


manifoldEdgeFlag = ...
    edgeCount == 2;


nonManifoldEdgeFlag = ...
    edgeCount > 2;


numberBoundaryEdges = ...
    sum(boundaryEdgeFlag);


numberNonManifoldEdges = ...
    sum(nonManifoldEdgeFlag);


if failOnNonManifold && ...
        numberNonManifoldEdges > 0

    error( ...
        ['fixNormalDirection: 检测到 %d 条 non-manifold 边，', ...
         '无法安全统一法向。'], ...
         numberNonManifoldEdges);

end


%% ========================================================================
% 9. 建立“面 - 面”方向约束
% ========================================================================
%
% 定义 orientationSign：
%
%       +1 = 保持原三角面方向
%
%       -1 = 翻转三角面方向
%
%
% 对公共边：
%
%       d1*s1 = -(d2*s2)
%
% 因此：
%
%       s2 = -d1*d2*s1
%

neighborFace = ...
    cell(numberFaces,1);


neighborRelation = ...
    cell(numberFaces,1);


for edgeIndex = 1:numberEdges

    if ~manifoldEdgeFlag(edgeIndex)

        continue;

    end


    occurrence = ...
        edgeOccurrence{edgeIndex};


    occurrence1 = ...
        occurrence(1);


    occurrence2 = ...
        occurrence(2);


    face1 = ...
        faceID(occurrence1);


    face2 = ...
        faceID(occurrence2);


    direction1 = ...
        directionSign(occurrence1);


    direction2 = ...
        directionSign(occurrence2);


    relation = ...
        -direction1 ...
        * ...
        direction2;


    neighborFace{face1}(end+1) = ...
        face2;


    neighborRelation{face1}(end+1) = ...
        relation;


    neighborFace{face2}(end+1) = ...
        face1;


    neighborRelation{face2}(end+1) = ...
        relation;

end


%% ========================================================================
% 10. BFS 遍历每一个连通区域
% ========================================================================

orientationSign = ...
    zeros(numberFaces,1);


componentID = ...
    zeros(numberFaces,1);


componentNumber = 0;


orientationConflictNumber = 0;


queue = ...
    zeros(numberFaces,1);


for startFace = 1:numberFaces

    if orientationSign(startFace) ~= 0

        continue;

    end


    componentNumber = ...
        componentNumber + 1;


    orientationSign(startFace) = 1;

    componentID(startFace) = ...
        componentNumber;


    queueStart = 1;

    queueEnd = 1;


    queue(queueEnd) = ...
        startFace;


    while queueStart <= queueEnd

        currentFace = ...
            queue(queueStart);


        queueStart = ...
            queueStart + 1;


        neighbours = ...
            neighborFace{currentFace};


        relations = ...
            neighborRelation{currentFace};


        for k = 1:numel(neighbours)

            nextFace = ...
                neighbours(k);


            expectedSign = ...
                relations(k) ...
                * ...
                orientationSign(currentFace);


            if orientationSign(nextFace) == 0

                orientationSign(nextFace) = ...
                    expectedSign;


                componentID(nextFace) = ...
                    componentNumber;


                queueEnd = ...
                    queueEnd + 1;


                queue(queueEnd) = ...
                    nextFace;


            elseif orientationSign(nextFace) ~= ...
                    expectedSign


                orientationConflictNumber = ...
                    orientationConflictNumber + 1;

            end

        end

    end

end


%% ========================================================================
% 11. 根据局部一致性结果翻转三角形
% ========================================================================

flipMask = ...
    false(numberFaces,1);


localFlip = ...
    orientationSign < 0;


flipMask(localFlip) = ...
    ~flipMask(localFlip);


F(localFlip,[2 3]) = ...
    F(localFlip,[3 2]);


numberLocalFlips = ...
    sum(localFlip);


%% ========================================================================
% 12. 判断每一个连通区域是否为开放区域
% ========================================================================

componentOpen = ...
    false(componentNumber,1);


for edgeIndex = 1:numberEdges

    if manifoldEdgeFlag(edgeIndex)

        continue;

    end


    occurrence = ...
        edgeOccurrence{edgeIndex};


    affectedFaces = ...
        unique( ...
        faceID(occurrence));


    affectedComponents = ...
        unique( ...
        componentID(affectedFaces));


    affectedComponents = ...
        affectedComponents( ...
        affectedComponents > 0);


    componentOpen(affectedComponents) = ...
        true;

end


%% ========================================================================
% 13. 根据整体有向体积判断外法向
% ========================================================================
%
% 对封闭三角表面：
%
%           1
% V = --------- Σ P1 · (P2 × P3)
%           6
%
%
% 对正确外法向：
%
%       V > 0
%
%
% 对整体反向：
%
%       V < 0
%

modelRange = ...
    max(V,[],1) ...
    - ...
    min(V,[],1);


referenceLength = ...
    max(modelRange);


if referenceLength <= 0

    referenceLength = 1;

end


volumeTolerance = ...
    1e-12 ...
    * ...
    referenceLength^3;


componentSignedVolume = ...
    zeros(componentNumber,1);


componentGlobalFlipped = ...
    false(componentNumber,1);


for component = 1:componentNumber

    faces = ...
        find( ...
        componentID == component);


    Fc = ...
        F(faces,:);


    P1 = ...
        V(Fc(:,1),:);


    P2 = ...
        V(Fc(:,2),:);


    P3 = ...
        V(Fc(:,3),:);


    % Evaluate volume near the component, avoiding cancellation after translation.
    volumeOrigin = mean(V(unique(Fc(:)),:),1);
    signedVolume = ...
        sum( ...
        sum( ...
        (P1-volumeOrigin) .* ...
        cross(P2-volumeOrigin,P3-volumeOrigin,2), ...
        2)) ...
        / ...
        6;


    componentSignedVolume(component) = ...
        signedVolume;


    %% ====================================================================
    % 13A. 封闭区域：优先使用有向体积
    %% ====================================================================

    if ~componentOpen(component) ...
            && ...
            abs(signedVolume) > volumeTolerance


        if signedVolume < 0

            F(faces,[2 3]) = ...
                F(faces,[3 2]);


            flipMask(faces) = ...
                ~flipMask(faces);


            componentGlobalFlipped(component) = ...
                true;


            componentSignedVolume(component) = ...
                -signedVolume;

        end


    %% ====================================================================
    % 13B. 开放区域或体积接近零：
    % 使用区域中心法作为备用
    %% ====================================================================

    else

        componentVertices = ...
            unique(Fc(:));


        regionCenter = ...
            mean( ...
            V(componentVertices,:), ...
            1);


        rawNormal = ...
            cross( ...
            P2-P1, ...
            P3-P1, ...
            2);


        normalMagnitude = ...
            sqrt( ...
            sum(rawNormal.^2,2));


        faceArea = ...
            0.5 ...
            * ...
            normalMagnitude;


        faceNormal = ...
            rawNormal ...
            ./ ...
            max( ...
            normalMagnitude, ...
            eps);


        faceCenter = ...
            (P1+P2+P3)/3;


        radialVector = ...
            faceCenter ...
            - ...
            regionCenter;


        outwardScore = ...
            sum( ...
            faceArea ...
            .* ...
            sum( ...
            faceNormal ...
            .* ...
            radialVector, ...
            2));


        if outwardScore < 0

            F(faces,[2 3]) = ...
                F(faces,[3 2]);


            flipMask(faces) = ...
                ~flipMask(faces);


            componentGlobalFlipped(component) = ...
                true;

        end

    end

end


%% ========================================================================
% 14. 将修正后的 F 写回 mesh
% ========================================================================

mesh.F = F;


%% ========================================================================
% 15. 重新计算面元几何
% ========================================================================

mesh = ...
    recalculateFaceGeometry(mesh);


%% ========================================================================
% 16. 保存法向检查信息
% ========================================================================

mesh.normalInfo.totalFaces = ...
    numberFaces;


mesh.normalInfo.connectedComponents = ...
    componentNumber;


mesh.normalInfo.boundaryEdges = ...
    numberBoundaryEdges;


mesh.normalInfo.nonManifoldEdges = ...
    numberNonManifoldEdges;


mesh.normalInfo.orientationConflicts = ...
    orientationConflictNumber;


mesh.normalInfo.localFlippedFaces = ...
    numberLocalFlips;


mesh.normalInfo.globallyFlippedComponents = ...
    sum(componentGlobalFlipped);


mesh.normalInfo.totalChangedFaces = ...
    sum(flipMask);


mesh.normalInfo.componentSignedVolume = ...
    componentSignedVolume;


mesh.normalInfo.isClosedSurface = ...
    numberBoundaryEdges == 0 ...
    && ...
    numberNonManifoldEdges == 0;


mesh.normalInfo.statusOK = ...
    numberNonManifoldEdges == 0 ...
    && ...
    orientationConflictNumber == 0;


%% ========================================================================
% 17. 输出报告
% ========================================================================

if verbose

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' Surface Normal Orientation Report\n');
    fprintf('============================================================\n');


    fprintf('Total faces              : %d\n', ...
        numberFaces);


    fprintf('Connected components     : %d\n', ...
        componentNumber);


    fprintf('\n');


    fprintf('Boundary edges           : %d\n', ...
        numberBoundaryEdges);


    fprintf('Non-manifold edges       : %d\n', ...
        numberNonManifoldEdges);


    fprintf('Orientation conflicts    : %d\n', ...
        orientationConflictNumber);


    fprintf('\n');


    fprintf('Local face flips         : %d\n', ...
        numberLocalFlips);


    fprintf('Global component flips   : %d\n', ...
        sum(componentGlobalFlipped));


    fprintf('Final changed faces      : %d\n', ...
        sum(flipMask));


    fprintf('\n');


    if mesh.normalInfo.isClosedSurface

        fprintf('Surface type             : CLOSED\n');

    else

        fprintf('Surface type             : OPEN / NON-MANIFOLD\n');

    end


    if mesh.normalInfo.statusOK

        fprintf('Normal orientation       : OK\n');

    else

        fprintf('Normal orientation       : WARNING\n');

    end


    fprintf('============================================================\n');

end

end



%% ========================================================================
% LOCAL FUNCTION
% 重新计算面元几何
%% ========================================================================
function mesh = recalculateFaceGeometry(mesh)

V = ...
    double(mesh.V);


F = ...
    double(mesh.F);


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
    sum(rawNormal.^2,2));


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
    (P1+P2+P3)/3;


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


modelLength = ...
    max(mesh.size);


if modelLength <= 0

    modelLength = 1;

end


areaTolerance = ...
    1e-14 ...
    * ...
    modelLength^2;


mesh.badFace = ...
    mesh.area <= ...
    areaTolerance;


mesh.numberBadFaces = ...
    sum(mesh.badFace);

end



%% ========================================================================
% LOCAL FUNCTION
% 安全读取 cfg 参数
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
