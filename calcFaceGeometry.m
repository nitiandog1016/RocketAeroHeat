function mesh = calcFaceGeometry(mesh)
%% =========================================================
%  calcFaceGeometry.m
%
%  功能：
%  根据 STL 的节点和三角面连接关系，
%  计算每一个三角面元的几何参数。
%
%
%  对第 i 个三角面：
%
%        P3
%        /\
%       /  \
%      /    \
%     /      \
%    P1------P2
%
%
%  定义：
%
%  e1 = P2 - P1
%  e2 = P3 - P1
%
%
%  三角形法向：
%
%        e1 × e2
%
%
%  三角形面积：
%
%                   1
%        A = ----------- |e1 × e2|
%                   2
%
%
%  单位法向：
%
%              e1 × e2
%        n = --------------
%             |e1 × e2|
%
%
%  面心：
%
%        C = (P1 + P2 + P3)/3
%
%
%  输入：
%
%  mesh.V
%      节点坐标
%
%  mesh.F
%      三角面连接关系
%
%
%  输出增加：
%
%  mesh.center
%      每个面元的面心坐标
%
%  mesh.normal
%      每个面元的单位法向
%
%  mesh.area
%      每个三角面元面积
%
%  mesh.totalArea
%      STL 总表面积
%
%  mesh.boundingBox
%      模型包围盒
%
%  mesh.size
%      X/Y/Z 三个方向的模型尺寸
%% =========================================================


%% ===================== 1. 提取节点和三角面 =====================

V = mesh.V;

F = mesh.F;


%% ===================== 2. 找到每个面元的三个顶点 =====================
% F(:,1)、F(:,2)、F(:,3)
% 分别保存每个三角形三个节点的编号。
%
% 根据这些编号直接从 V 中提取三维坐标。

P1 = V(F(:,1), :);

P2 = V(F(:,2), :);

P3 = V(F(:,3), :);


%% ===================== 3. 计算两条三角形边向量 =====================

e1 = P2 - P1;

e2 = P3 - P1;


%% ===================== 4. 计算叉乘 =====================
% cross(...,2)
%
% 表示沿矩阵的第 2 个维度，
% 对每一行分别进行三维叉乘。
%
% 得到的 normalRaw：
%
% 方向：
% 三角面法向方向
%
% 模长：
% |e1 × e2| = 2A

normalRaw = cross(e1, e2, 2);


%% ===================== 5. 计算叉乘模长 =====================

normalMagnitude = sqrt(sum(normalRaw.^2, 2));


%% ===================== 6. 检查退化三角形 =====================
% 如果三个点共线，或者两个点完全重合：
%
% |e1 × e2| = 0
%
% 那么三角形面积为 0，
% 法向也无法定义。
%
% 这类三角形属于几何缺陷。

tolerance = 1e-14;

badFace = normalMagnitude < tolerance;

numberBadFaces = sum(badFace);


if numberBadFaces > 0

    warning(['Detected %d degenerate triangular faces ', ...
             'with nearly zero area.'], ...
             numberBadFaces);

end


%% ===================== 7. 计算三角形面积 =====================

area = 0.5 * normalMagnitude;


%% ===================== 8. 计算单位法向 =====================
% 为避免除以 0，
% 先将异常面元对应模长临时设置为 1。
%
% 异常面元最终法向设为 NaN。

safeMagnitude = normalMagnitude;

safeMagnitude(badFace) = 1;


normal = normalRaw ./ safeMagnitude;


% 将退化面元的法向设为 NaN
normal(badFace,:) = NaN;


%% ===================== 9. 计算面心 =====================
% 三角形几何中心：
%
% C = (P1 + P2 + P3)/3

center = (P1 + P2 + P3) / 3;


%% ===================== 10. 计算总表面积 =====================

totalArea = sum(area);


%% ===================== 11. 计算模型包围盒 =====================
% boundingBox：
%
% 第一行：Xmin Xmax
% 第二行：Ymin Ymax
% 第三行：Zmin Zmax

xmin = min(V(:,1));
xmax = max(V(:,1));

ymin = min(V(:,2));
ymax = max(V(:,2));

zmin = min(V(:,3));
zmax = max(V(:,3));


boundingBox = [
    xmin xmax;
    ymin ymax;
    zmin zmax
];


%% ===================== 12. 计算模型三个方向尺寸 =====================

modelSize = [
    xmax - xmin, ...
    ymax - ymin, ...
    zmax - zmin
];


%% ===================== 13. 保存计算结果 =====================

mesh.center = center;

mesh.normal = normal;

mesh.area = area;

mesh.badFace = badFace;

mesh.numberBadFaces = numberBadFaces;

mesh.totalArea = totalArea;

mesh.boundingBox = boundingBox;

mesh.size = modelSize;


%% ===================== 14. 输出计算结果 =====================

fprintf('Face centers calculated.\n');

fprintf('Face normals calculated.\n');

fprintf('Face areas calculated.\n');

fprintf('Degenerate faces      : %d\n', ...
    numberBadFaces);

end