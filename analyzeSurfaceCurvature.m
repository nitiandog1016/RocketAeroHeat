function mesh = analyzeSurfaceCurvature(mesh, cfg)
%% =========================================================
% analyzeSurfaceCurvature.m
%
% 功能：
%
% 对 STL 三角面网格进行离散曲率分析。
%
% 核心思想：
%
% 如果两个相邻三角面的法向几乎相同：
%
%       n1 ≈ n2
%
% 说明局部表面比较平缓，曲率较小，
% 可以使用比较大的三角面。
%
%
% 如果两个相邻三角面的法向变化明显：
%
%       n1 != n2
%
% 说明局部表面弯曲明显，
% 需要使用更小的三角面。
%
%
% 离散曲率近似：
%
%                 Δtheta
%       kappa ≈ ----------
%                  Δs
%
%
% 其中：
%
% Δtheta：
% 两个相邻面法向夹角 [rad]
%
% Δs：
% 两面公共边长度 [m]
%
%
% 输入：
%
% mesh.V
%     节点坐标
%
% mesh.F
%     三角面连接关系
%
% mesh.normal
%     每个三角面的单位法向
%
%
% 输出增加：
%
% mesh.curvatureMean
%     每个面元平均离散曲率 [1/m]
%
% mesh.curvatureMax
%     每个面元最大离散曲率 [1/m]
%
% mesh.targetSize
%     根据曲率计算得到的建议面元尺寸 [m]
%
% mesh.featureFace
%     是否位于明显几何特征附近
%
% mesh.featureEdges
%     检测到的特征边节点编号
%
%% =========================================================


%% ===================== 1. 读取网格 =====================

V = mesh.V;
F = mesh.F;
N = mesh.normal;

numberFaces = size(F,1);


%% ===================== 2. 建立所有三角形边 =====================
%
% 每个三角形包含三条边：
%
% edge 1 : node1 - node2
% edge 2 : node2 - node3
% edge 3 : node3 - node1
%
% 把所有三角形的边全部列出来。

edges = [
    F(:,[1 2]);
    F(:,[2 3]);
    F(:,[3 1])
];


%% ===================== 3. 保存每条边属于哪个面元 =====================
%
% 第一批边对应面元 1...N
% 第二批边也对应面元 1...N
% 第三批边同样如此。

faceID = repmat((1:numberFaces)', 3, 1);


%% ===================== 4. 对边的节点编号排序 =====================
%
% 例如：
%
% [5 10]
%
% 和：
%
% [10 5]
%
% 实际上是同一条边。
%
% 因此把每条边统一变成：
%
% [minNode maxNode]

edgesSorted = sort(edges, 2);


%% ===================== 5. 寻找重复边 =====================
%
% 对封闭三角表面：
%
% 内部公共边应该同时属于两个三角面。
%
% 如果一条边只属于一个面：
% 说明它可能是自由边 / 开口边界。
%
% 如果属于超过两个面：
% 可能存在 non-manifold 几何问题。

[uniqueEdges, ~, edgeGroup] = unique(edgesSorted, 'rows');


%% ===================== 6. 统计每条边出现次数 =====================

edgeCount = accumarray(edgeGroup, 1);


%% ===================== 7. 找到每条边相邻的面元编号 =====================
%
% 对正常公共边：
%
% edgeCount = 2
%
% 于是可以得到：
%
% face1
% face2

firstFace = accumarray( ...
    edgeGroup, ...
    faceID, ...
    [], ...
    @min);

secondFace = accumarray( ...
    edgeGroup, ...
    faceID, ...
    [], ...
    @max);


%% ===================== 8. 提取正常公共边 =====================

internalEdgeFlag = ...
    (edgeCount == 2) & ...
    (firstFace ~= secondFace);

internalEdges = uniqueEdges(internalEdgeFlag,:);

face1 = firstFace(internalEdgeFlag);

face2 = secondFace(internalEdgeFlag);


%% ===================== 9. 计算公共边长度 =====================

P1 = V(internalEdges(:,1),:);

P2 = V(internalEdges(:,2),:);

edgeVector = P2 - P1;

edgeLength = sqrt(sum(edgeVector.^2,2));


%% ===================== 10. 获取公共边两侧三角面法向 =====================

N1 = N(face1,:);

N2 = N(face2,:);


%% ===================== 11. 检查有效数据 =====================
%
% 如果某些三角面是退化三角形，
% 法向可能是 NaN。
%
% 这些边暂时不参与曲率计算。

validEdge = ...
    all(isfinite(N1),2) & ...
    all(isfinite(N2),2) & ...
    edgeLength > eps;


face1Valid = face1(validEdge);

face2Valid = face2(validEdge);

edgeLengthValid = edgeLength(validEdge);

N1Valid = N1(validEdge,:);

N2Valid = N2(validEdge,:);

internalEdgesValid = internalEdges(validEdge,:);


%% ===================== 12. 计算相邻面法向夹角 =====================
%
% dot(n1,n2) = cos(theta)
%
% 根据共享边的绕序修正局部法向符号，保留大于 90 度的真实折角。

dotNormal = sum(N1Valid .* N2Valid, 2);

direction1 = sharedEdgeDirection(F(face1Valid,:),internalEdgesValid);
direction2 = sharedEdgeDirection(F(face2Valid,:),internalEdgesValid);
dotNormal = dotNormal .* (-direction1 .* direction2);


% 数值误差可能让 dot 略微超过 1，
% 所以强制限制到 [-1,1]。

dotNormal = min(1, max(-1, dotNormal));


edgeAngle = acos(dotNormal);      % rad

edgeAngleDeg = edgeAngle * 180/pi;


%% ===================== 13. 计算离散曲率 =====================
%
% kappa ≈ Δtheta / Δs
%
% 单位：
%
% rad/m ≈ 1/m

edgeCurvature = ...
    edgeAngle ./ edgeLengthValid;


%% ===================== 14. 将边曲率映射回三角面 =====================
%
% 每个三角面通常有三条边。
%
% 这里分别计算：
%
% 平均曲率
% 最大曲率

allAdjacentFaces = [
    face1Valid;
    face2Valid
];

allCurvature = [
    edgeCurvature;
    edgeCurvature
];


%% 平均曲率

curvatureMean = accumarray( ...
    allAdjacentFaces, ...
    allCurvature, ...
    [numberFaces 1], ...
    @mean, ...
    0);


%% 最大曲率

curvatureMax = accumarray( ...
    allAdjacentFaces, ...
    allCurvature, ...
    [numberFaces 1], ...
    @max, ...
    0);


%% ===================== 15. 判断几何特征边 =====================
%
% 如果两个相邻三角面法向夹角超过：
%
% cfg.remesh.featureAngle
%
% 则认为这里存在明显折角或特征。

featureEdgeFlag = ...
    edgeAngleDeg >= cfg.remesh.featureAngle;


featureEdges = ...
    internalEdgesValid(featureEdgeFlag,:);


%% ===================== 16. 找出特征边附近的面 =====================

featureFace = false(numberFaces,1);

featureFace(face1Valid(featureEdgeFlag)) = true;

featureFace(face2Valid(featureEdgeFlag)) = true;


%% ===================== 17. 根据曲率计算目标面元尺寸 =====================
%
% 理想情况：
%
% 一个三角面所跨越的法向变化不要超过：
%
% cfg.remesh.curvatureAngleTarget
%
%
% 如果：
%
% kappa ≈ Δtheta / h
%
% 则：
%
% h ≈ Δtheta_target / kappa
%
%
% 曲率越大：
%
% kappa ↑
%
% 则：
%
% h ↓
%
% 即高曲率区域自动加密。

thetaTarget = ...
    cfg.remesh.curvatureAngleTarget ...
    * pi / 180;


sensitivity = ...
    cfg.remesh.curvatureSensitivity;


targetSize = ...
    cfg.remesh.hMax ...
    * ones(numberFaces,1);


curvedFace = curvatureMax > 1e-12;


targetSize(curvedFace) = ...
    thetaTarget ./ ...
    (sensitivity .* curvatureMax(curvedFace));


%% ===================== 18. 限制最大最小尺寸 =====================

targetSize = max( ...
    targetSize, ...
    cfg.remesh.hMin);


targetSize = min( ...
    targetSize, ...
    cfg.remesh.hMax);


%% ===================== 19. 特征区域进一步加密 =====================
%
% 特征边附近使用：
%
% featureSizeFactor × hMin

featureTargetSize = ...
    cfg.remesh.featureSizeFactor ...
    * cfg.remesh.hMin;


featureTargetSize = min( ...
    featureTargetSize, ...
    cfg.remesh.hMax);


targetSize(featureFace) = min( ...
    targetSize(featureFace), ...
    featureTargetSize);


%% ===================== 20. 保存结果 =====================

mesh.curvatureMean = curvatureMean;

mesh.curvatureMax = curvatureMax;

mesh.targetSize = targetSize;

mesh.featureFace = featureFace;

mesh.featureEdges = featureEdges;


%% ===================== 21. 保存边分析信息 =====================

mesh.edgeInfo.edges = internalEdgesValid;

mesh.edgeInfo.angleDeg = edgeAngleDeg;

mesh.edgeInfo.curvature = edgeCurvature;


%% ===================== 22. 输出信息 =====================

fprintf('\nSurface curvature analysis completed.\n');

fprintf('Internal shared edges     : %d\n', ...
    size(internalEdgesValid,1));

fprintf('Feature edges             : %d\n', ...
    size(featureEdges,1));

fprintf('Feature faces             : %d\n', ...
    sum(featureFace));

fprintf('Maximum discrete curvature: %.6e 1/m\n', ...
    max(curvatureMax));

fprintf('Minimum target size       : %.6e m\n', ...
    min(targetSize));

fprintf('Maximum target size       : %.6e m\n', ...
    max(targetSize));

end

function direction = sharedEdgeDirection(faces,edges)
forward = false(size(faces,1),1);
for k = 1:3
    j = mod(k,3)+1;
    forward = forward | ...
        (faces(:,k) == edges(:,1) & faces(:,j) == edges(:,2));
end
direction = 2*double(forward)-1;
end
