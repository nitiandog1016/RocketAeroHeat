function plotMeshQuality(mesh, quality, cfg)
%% =========================================================
% plotMeshQuality.m
%
% 功能：
%
% 将三角网格质量直接映射到 STL 表面。
%
%
% 左图：
%
% Aspect Ratio
%
% AR ≈ 1
%     三角形质量很好
%
% AR 越大
%     三角形越细长
%
%
% 右图：
%
% Minimum Triangle Angle
%
% 接近 60°
%     最理想
%
% 很小，例如：
%
% < 10°
%
%     表明存在非常尖锐的瘦三角形。
%
%% =========================================================


%% ===================== 1. 判断是否绘图 =====================

if isfield(cfg.plot,'showMeshQuality')

    if ~cfg.plot.showMeshQuality
        return;
    end

end


%% ===================== 2. 创建图窗 =====================

figure( ...
    'Color','w', ...
    'Name','Triangle Mesh Quality');


%% =========================================================
% 左侧：Aspect Ratio
%% =========================================================

subplot(1,2,1);


patch( ...
    'Faces', mesh.F, ...
    'Vertices', mesh.V, ...
    'FaceVertexCData', quality.aspectRatio, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none');


axis equal;
axis tight;

box on;
grid on;

view(35,25);


xlabel('X / m');
ylabel('Y / m');
zlabel('Z / m');


title('Triangle Aspect Ratio');


cb1 = colorbar;

ylabel(cb1,'Aspect Ratio');


%% ===================== Aspect Ratio 显示范围 =====================
%
% 如果某几个面 AR 极大，
% 直接按最大值绘图会让其他区域颜色全部挤在一起。
%
% 因此云图上限设置为：
%
% 2 × 允许最大 AR

ARmaxPlot = ...
    max( ...
    2, ...
    2 * cfg.remesh.maxAspectRatio);


caxis([1 ARmaxPlot]);


colormap(gca,parula);


%% ===================== 字体 =====================

ax1 = gca;

if isfield(cfg.plot,'fontName')

    ax1.FontName = cfg.plot.fontName;

end

if isfield(cfg.plot,'fontSize')

    ax1.FontSize = cfg.plot.fontSize;

end


%% =========================================================
% 右侧：最小内角
%% =========================================================

subplot(1,2,2);


patch( ...
    'Faces', mesh.F, ...
    'Vertices', mesh.V, ...
    'FaceVertexCData', quality.minAngle, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none');


axis equal;
axis tight;

box on;
grid on;

view(35,25);


xlabel('X / m');
ylabel('Y / m');
zlabel('Z / m');


title('Minimum Triangle Angle');


cb2 = colorbar;

ylabel(cb2,'Minimum Angle / deg');


%% ===================== 最小角显示范围 =====================
%
% 理想等边三角形：
%
% 60°
%
% 因此显示：
%
% 0 ~ 60°

caxis([0 60]);


colormap(gca,parula);


%% ===================== 字体 =====================

ax2 = gca;

if isfield(cfg.plot,'fontName')

    ax2.FontName = cfg.plot.fontName;

end

if isfield(cfg.plot,'fontSize')

    ax2.FontSize = cfg.plot.fontSize;

end


%% ===================== 总标题 =====================

sgtitle( ...
    sprintf( ...
    'Mesh Quality: %d Faces, %.2f %% Poor Faces', ...
    quality.numberFaces, ...
    quality.badFacePercent));


end