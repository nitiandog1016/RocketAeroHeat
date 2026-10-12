function quality = checkMeshQuality(mesh, cfg)
%% =========================================================
% checkMeshQuality.m
%
% 功能：
%
% 对三角表面网格进行质量评价。
%
% 主要检查：
%
% 1. 三角形面积
% 2. 三条边长度
% 3. 三个内角
% 4. 最小内角
% 5. 长宽比 Aspect Ratio
% 6. 最长边 / 最短边
% 7. 差质量三角面数量
%
%
% 理想三角形：
%
%        /\
%       /  \
%      /    \
%     /______\
%
% 三个角接近：
%
% 60°
%
% Aspect Ratio 接近：
%
% 1
%
%
% 面元越细长：
%
% AR 越大
%
% 最小角越小。
%
%% =========================================================


%% ===================== 1. 获取节点和三角形 =====================

V = mesh.V;

F = mesh.F;

numberFaces = size(F,1);


%% ===================== 2. 提取三个顶点 =====================

P1 = V(F(:,1),:);

P2 = V(F(:,2),:);

P3 = V(F(:,3),:);


%% ===================== 3. 计算三条边长度 =====================
%
% 定义：
%
% a = |P2-P3|
% b = |P1-P3|
% c = |P1-P2|
%
% 其中：
%
% a 对应 P1 对面的边
% b 对应 P2 对面的边
% c 对应 P3 对面的边

a = sqrt(sum((P2-P3).^2,2));

b = sqrt(sum((P1-P3).^2,2));

c = sqrt(sum((P1-P2).^2,2));


%% ===================== 4. 最大边和最小边 =====================

edgeMatrix = [a b c];

maxEdge = max(edgeMatrix,[],2);

minEdge = min(edgeMatrix,[],2);


%% ===================== 5. 最长边 / 最短边 =====================

edgeRatio = ...
    maxEdge ./ max(minEdge,eps);


%% ===================== 6. 三角形面积 =====================

if isfield(mesh,'area')

    area = mesh.area;

else

    normalRaw = cross(P2-P1,P3-P1,2);

    area = 0.5 * ...
        sqrt(sum(normalRaw.^2,2));

end


%% ===================== 7. 利用余弦定理计算三个内角 =====================
%
% angle1 位于 P1
%
% cos(angle1) =
%
%     b^2 + c^2 - a^2
% -------------------------
%          2bc

cosAngle1 = ...
    (b.^2 + c.^2 - a.^2) ./ ...
    max(2 .* b .* c, eps);


cosAngle2 = ...
    (a.^2 + c.^2 - b.^2) ./ ...
    max(2 .* a .* c, eps);


cosAngle3 = ...
    (a.^2 + b.^2 - c.^2) ./ ...
    max(2 .* a .* b, eps);


%% ===================== 8. 防止浮点误差 =====================

cosAngle1 = min(1,max(-1,cosAngle1));

cosAngle2 = min(1,max(-1,cosAngle2));

cosAngle3 = min(1,max(-1,cosAngle3));


%% ===================== 9. 转换成角度 =====================

angle1 = acos(cosAngle1) * 180/pi;

angle2 = acos(cosAngle2) * 180/pi;

angle3 = acos(cosAngle3) * 180/pi;


angleMatrix = [
    angle1 ...
    angle2 ...
    angle3
];


%% ===================== 10. 最小/最大角 =====================

minAngle = min(angleMatrix,[],2);

maxAngle = max(angleMatrix,[],2);


%% ===================== 11. 计算三角形内切圆半径 =====================
%
% 三角形内切圆半径：
%
%             2A
%       r = -------
%            a+b+c

perimeter = a + b + c;

inRadius = ...
    2 .* area ./ ...
    max(perimeter,eps);


%% ===================== 12. Aspect Ratio =====================
%
% 这里采用：
%
%                    Lmax
% AR = --------------------------------
%       2 * sqrt(3) * r_in
%
%
% 对于理想等边三角形：
%
% AR = 1
%
% AR 越大：
%
% 三角形越细长。

aspectRatio = ...
    maxEdge ./ ...
    max(2*sqrt(3).*inRadius,eps);


%% ===================== 13. 判断退化面 =====================

degenerateFace = ...
    (area < 1e-14) | ...
    ~isfinite(aspectRatio);


%% ===================== 14. 判断长宽比差的面元 =====================

badAspectRatio = ...
    aspectRatio > ...
    cfg.remesh.maxAspectRatio;


%% ===================== 15. 判断最小角太小的面元 =====================

badAngle = ...
    minAngle < ...
    cfg.remesh.minAngle;


%% ===================== 16. 加入原来的坏面信息 =====================

badOriginal = false(numberFaces,1);

if isfield(mesh,'badFace')

    badOriginal = mesh.badFace;

end


%% ===================== 17. 综合差质量判断 =====================

badFace = ...
    degenerateFace | ...
    badAspectRatio | ...
    badAngle | ...
    badOriginal;


%% ===================== 18. 保存逐面质量结果 =====================

quality.area = area;

quality.edgeLength = edgeMatrix;

quality.maxEdge = maxEdge;

quality.minEdge = minEdge;

quality.edgeRatio = edgeRatio;

quality.angle = angleMatrix;

quality.minAngle = minAngle;

quality.maxAngle = maxAngle;

quality.aspectRatio = aspectRatio;

quality.degenerateFace = degenerateFace;

quality.badAspectRatio = badAspectRatio;

quality.badAngle = badAngle;

quality.badFace = badFace;


%% ===================== 19. 保存统计信息 =====================

quality.numberFaces = numberFaces;

quality.numberBadFaces = sum(badFace);

quality.numberBadAR = sum(badAspectRatio);

quality.numberBadAngle = sum(badAngle);

quality.numberDegenerate = sum(degenerateFace);


validArea = area(area > 0 & isfinite(area));


if ~isempty(validArea)

    quality.minArea = min(validArea);

    quality.maxArea = max(validArea);

    quality.meanArea = mean(validArea);

else

    quality.minArea = NaN;

    quality.maxArea = NaN;

    quality.meanArea = NaN;

end


validAR = aspectRatio(isfinite(aspectRatio));


if ~isempty(validAR)

    quality.maxAspectRatio = max(validAR);

    quality.meanAspectRatio = mean(validAR);

else

    quality.maxAspectRatio = Inf;

    quality.meanAspectRatio = Inf;

end


validAngle = minAngle(isfinite(minAngle));


if ~isempty(validAngle)

    quality.minimumAngle = min(validAngle);

    quality.meanMinimumAngle = mean(validAngle);

else

    quality.minimumAngle = NaN;

    quality.meanMinimumAngle = NaN;

end


%% ===================== 20. 差质量比例 =====================

quality.badFacePercent = ...
    100 * quality.numberBadFaces ...
    / max(numberFaces,1);


%% ===================== 21. 输出统计信息 =====================

fprintf('\n');
fprintf('=============================================\n');
fprintf(' Triangle Mesh Quality Report\n');
fprintf('=============================================\n');

fprintf('Total faces               : %d\n', ...
    numberFaces);

fprintf('\n');

fprintf('Minimum area              : %.6e m^2\n', ...
    quality.minArea);

fprintf('Maximum area              : %.6e m^2\n', ...
    quality.maxArea);

fprintf('Mean area                 : %.6e m^2\n', ...
    quality.meanArea);

fprintf('\n');

fprintf('Minimum triangle angle    : %.3f deg\n', ...
    quality.minimumAngle);

fprintf('Mean minimum angle        : %.3f deg\n', ...
    quality.meanMinimumAngle);

fprintf('\n');

fprintf('Maximum aspect ratio      : %.3f\n', ...
    quality.maxAspectRatio);

fprintf('Mean aspect ratio         : %.3f\n', ...
    quality.meanAspectRatio);

fprintf('\n');

fprintf('AR > %.2f                 : %d\n', ...
    cfg.remesh.maxAspectRatio, ...
    quality.numberBadAR);

fprintf('Angle < %.2f deg          : %d\n', ...
    cfg.remesh.minAngle, ...
    quality.numberBadAngle);

fprintf('Degenerate faces          : %d\n', ...
    quality.numberDegenerate);

fprintf('\n');

fprintf('Total poor-quality faces  : %d\n', ...
    quality.numberBadFaces);

fprintf('Poor-quality percentage   : %.3f %%\n', ...
    quality.badFacePercent);

fprintf('=============================================\n');


end