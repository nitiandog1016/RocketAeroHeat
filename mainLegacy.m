%% ========================================================================
% main.m
%
% Rocket AeroHeat
%
% 三维箭体工程气动 / 气动热计算程序
%
% ========================================================================
%
% 主程序流程：
%   1. 读取 STL，在统一窗口中确认箭体朝向并重建表面网格；
%   2. 统一外法向、检查网格，建立从鼻尖向尾部的表面坐标；
%   3. 在同一窗口设置飞行条件，计算面元压力与气动力；
%   4. 内置工程外流、局部激波、三维表面准流线、边界层及转捩；
%   5. 估算驻点、偏置球鼻和表面热流，在同一窗口显示算法与适用范围；
%   6. 用户填写材料参数后，可求固定或多站工况的壁温历程，
%      并查看累计热量、峰值热流、已求解面积和能量守恒核对。
% 壁温和累计热量由 GUI 按需计算，不在主程序启动时臆造材料参数。
%
%
% ========================================================================
%
% 当前内部坐标标准：
%
%       +X = 鼻尖 / 飞行前向
%
%       +Y = 横向
%
%       +Z = 向上
%
%
% ========================================================================
%
% 主程序结束以后重要变量：
%
%       cfg
%
%       rawMeshOriginal
%
%       rawMesh
%
%       mesh
%
%       quality
%
%       flight
%
%       atmosphere
%
%       freestream
%
%
% 其中：
%
%       mesh.param
%
% 保存半结构化截面拓扑。
%
%
%       mesh.surfaceCoord
%           x / r / theta / s 坐标。
%       mesh.heating
%           当前飞行工况的工程热流与算法分区。
%       mesh.wallThermal / mesh.aerothermalExposure
%           用户在 GUI 中完成材料和时间设置后生成。
%
%% ========================================================================



%% ========================================================================
% 0. MATLAB 初始化
% ========================================================================

clc;

clear;

close all;



%% ========================================================================
% 1. 工程根目录
% ========================================================================

projectRoot = ...
    fileparts( ...
        mfilename('fullpath'));


addpath( ...
    genpath(projectRoot));



%% ========================================================================
% 2. 启动信息
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('                    Rocket AeroHeat\n');
fprintf('          Aero-Thermal Engineering Code\n');
fprintf('============================================================\n');

fprintf('\n');



%% ========================================================================
% STEP 1 / 9
% 读取配置
% ========================================================================

fprintf('============================================================\n');
fprintf('Step 1 / 9 : Reading configuration\n');
fprintf('============================================================\n');


cfg = ...
    createCaseConfig();



%% ========================================================================
% 3. STL 文件检查
% ========================================================================

if ~isfield(cfg,'geometry') ...
        || ...
        ~isfield(cfg.geometry,'stlFile')

    error( ...
        ['cfg.geometry.stlFile 不存在。' newline ...
         '请检查 createCaseConfig.m。']);

end


if ~isfile(cfg.geometry.stlFile)

    error( ...
        ['找不到 STL 文件：' newline ...
         cfg.geometry.stlFile]);

end



%% ========================================================================
% STEP 2 / 9
% 读取原始 STL
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 2 / 9 : Loading original STL geometry\n');
fprintf('============================================================\n');


rawMeshOriginal = ...
    loadSTLMesh(cfg);



%% ========================================================================
% 4. 原始 STL 几何信息
% ========================================================================

rawMeshOriginal = ...
    calcFaceGeometry( ...
        rawMeshOriginal);



%% ========================================================================
% 5. 原始几何报告
% ========================================================================

fprintf('\n');

fprintf('------------------------------------------------------------\n');
fprintf(' Original STL Coordinate System\n');
fprintf('------------------------------------------------------------\n');


fprintf('Vertices                 : %d\n', ...
    size(rawMeshOriginal.V,1));


fprintf('Triangles                : %d\n', ...
    size(rawMeshOriginal.F,1));


fprintf('Surface area             : %.8e m^2\n', ...
    rawMeshOriginal.totalArea);


fprintf('\n');


fprintf('Original dimensions:\n');


fprintf('X                         : %.6f m\n', ...
    rawMeshOriginal.size(1));


fprintf('Y                         : %.6f m\n', ...
    rawMeshOriginal.size(2));


fprintf('Z                         : %.6f m\n', ...
    rawMeshOriginal.size(3));


fprintf('\n');


fprintf('Current STL orientation is NOT assumed correct.\n');


fprintf('\n');


fprintf('Required internal convention:\n');

fprintf('    Body axis             = X\n');

fprintf('    Nose direction        = +X\n');

fprintf('    Lateral direction     = Y\n');

fprintf('    Up direction          = Z\n');


fprintf('------------------------------------------------------------\n');



%% ========================================================================
% STEP 3 / 9
% 几何朝向 + Surface Mesh GUI
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 3 / 9 : Geometry orientation and surface meshing\n');
fprintf('============================================================\n');


% 创建唯一主窗口，后续页面均嵌入 appHost。
[appFigure,appHost] = createAeroHeatWindow();
setAeroHeatStage(appFigure,1);

[cfg, ...
 rawMesh, ...
 mesh, ...
 quality, ...
 meshAction] = ...
    meshGUI( ...
        cfg, ...
        rawMeshOriginal, appHost);



%% ========================================================================
% 6. 用户取消几何/网格前处理
% ========================================================================

if strcmpi( ...
        meshAction, ...
        'cancel')

    fprintf('\n');

    fprintf('============================================================\n');
    fprintf(' Geometry / mesh preprocessing cancelled.\n');
    fprintf(' Main program terminated.\n');
    fprintf('============================================================\n');

    fprintf('\n');


    if isgraphics(appFigure), delete(appFigure); end
    return;

end



%% ========================================================================
% 7. 几何朝向确认检查
% ========================================================================

if ~isfield(cfg.geometry,'isAligned') ...
        || ...
        ~cfg.geometry.isAligned

    error( ...
        ['几何朝向尚未确认。' newline ...
         '必须满足：' newline ...
         'Body axis = X' newline ...
         'Nose = +X']);

end



%% ========================================================================
% 8. 检查返回网格
% ========================================================================

if isempty(rawMesh)

    error( ...
        'meshGUI 没有返回对齐后的 rawMesh。');

end


if isempty(mesh)

    error( ...
        'meshGUI 没有返回最终计算网格 mesh。');

end



%% ========================================================================
% 9. 输出几何旋转信息
% ========================================================================

fprintf('\n');

fprintf('------------------------------------------------------------\n');
fprintf(' Geometry Alignment\n');
fprintf('------------------------------------------------------------\n');


fprintf('Alignment method         : %s\n', ...
    cfg.geometry.alignmentMethod);


fprintf('\n');


fprintf('Rotation matrix:\n');

disp( ...
    cfg.geometry.rotationMatrix);


fprintf('Rotation center:\n');

disp( ...
    cfg.geometry.rotationCenter);


fprintf('\n');


fprintf('Aligned dimensions:\n');


fprintf('X                         : %.6f m\n', ...
    rawMesh.size(1));


fprintf('Y                         : %.6f m\n', ...
    rawMesh.size(2));


fprintf('Z                         : %.6f m\n', ...
    rawMesh.size(3));


fprintf('\n');


fprintf('Internal body axis       : +X\n');

fprintf('Internal nose direction  : +X\n');


fprintf('------------------------------------------------------------\n');



%% ========================================================================
% STEP 4 / 9
% 表面外法向统一
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 4 / 9 : Surface normal orientation\n');
fprintf('============================================================\n');


mesh = ...
    fixNormalDirection( ...
        mesh, ...
        cfg);



%% ========================================================================
% STEP 5 / 9
% 最终表面网格验证
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 5 / 9 : Final surface mesh verification\n');
fprintf('============================================================\n');



%% ========================================================================
% 10. 确保最终面元几何完整
% ========================================================================

if ~isfield(mesh,'center') ...
        || ...
        ~isfield(mesh,'normal') ...
        || ...
        ~isfield(mesh,'area')

    mesh = ...
        calcFaceGeometry(mesh);

end



%% ========================================================================
% 11. 最终质量检查
% ========================================================================

quality = ...
    checkMeshQuality( ...
        mesh, ...
        cfg);



%% ========================================================================
% 12. 最终拓扑检查
% ========================================================================

if isfield(mesh,'normalInfo')

    boundaryEdges = ...
        mesh.normalInfo.boundaryEdges;


    nonManifoldEdges = ...
        mesh.normalInfo.nonManifoldEdges;


else

    [boundaryEdges, ...
     nonManifoldEdges] = ...
        localCheckTopology( ...
            mesh.F);

end



%% ========================================================================
% 13. 表面积变化检查
%
% 刚体旋转不能改变面积。
%
% 因此：
%
% rawMeshOriginal → rawMesh
%
% 面积变化应该接近机器精度。
%
%
% rawMesh → mesh
%
% 才是真正由重新面元化引起的几何误差。
%% ========================================================================

originalArea = ...
    rawMeshOriginal.totalArea;


alignedArea = ...
    rawMesh.totalArea;


finalArea = ...
    mesh.totalArea;


alignmentAreaChange = ...
    100 ...
    * ...
    (alignedArea-originalArea) ...
    / ...
    max(originalArea,eps);


remeshAreaChange = ...
    100 ...
    * ...
    (finalArea-alignedArea) ...
    / ...
    max(alignedArea,eps);



%% ========================================================================
% 14. 最终网格报告
% ========================================================================

fprintf('\n');

fprintf('------------------------------------------------------------\n');
fprintf(' Final Surface Mesh\n');
fprintf('------------------------------------------------------------\n');


fprintf('Original STL vertices    : %d\n', ...
    size(rawMeshOriginal.V,1));


fprintf('Final vertices           : %d\n', ...
    size(mesh.V,1));


fprintf('\n');


fprintf('Original STL triangles   : %d\n', ...
    size(rawMeshOriginal.F,1));


fprintf('Final triangles          : %d\n', ...
    size(mesh.F,1));


fprintf('\n');


fprintf('Original surface area    : %.8e m^2\n', ...
    originalArea);


fprintf('Aligned surface area     : %.8e m^2\n', ...
    alignedArea);


fprintf('Final mesh area          : %.8e m^2\n', ...
    finalArea);


fprintf('\n');


fprintf('Rotation area change     : %.8e %%\n', ...
    alignmentAreaChange);


fprintf('Remesh area change       : %.4f %%\n', ...
    remeshAreaChange);


fprintf('\n');


fprintf('Mean Aspect Ratio        : %.4f\n', ...
    quality.meanAspectRatio);


fprintf('Maximum Aspect Ratio     : %.4f\n', ...
    quality.maxAspectRatio);


fprintf('Mean minimum angle       : %.4f deg\n', ...
    quality.meanMinimumAngle);


fprintf('Minimum triangle angle   : %.4f deg\n', ...
    quality.minimumAngle);


fprintf('Poor-quality faces       : %d\n', ...
    quality.numberBadFaces);


fprintf('Poor-quality percentage  : %.4f %%\n', ...
    quality.badFacePercent);


fprintf('\n');


fprintf('Boundary edges           : %d\n', ...
    boundaryEdges);


fprintf('Non-manifold edges       : %d\n', ...
    nonManifoldEdges);



%% ========================================================================
% 15. 拓扑判断
% ========================================================================

if boundaryEdges == 0 ...
        && ...
        nonManifoldEdges == 0

    fprintf('Surface topology         : OK\n');


else

    fprintf('Surface topology         : WARNING\n');


    warning( ...
        ['最终表面网格存在拓扑问题：' newline ...
         'Boundary edges = %d' newline ...
         'Non-manifold edges = %d'], ...
         boundaryEdges, ...
         nonManifoldEdges);

end



%% ========================================================================
% 16. 几何保持判断
% ========================================================================

if isfinite(remeshAreaChange)

    if abs(remeshAreaChange) <= 1

        fprintf('Geometry preservation    : Very good\n');


    elseif abs(remeshAreaChange) <= 3

        fprintf('Geometry preservation    : Acceptable\n');


    else

        fprintf('Geometry preservation    : WARNING\n');


        warning( ...
            ['重新面元化前后表面积变化为 %.3f%%。' newline ...
             '建议检查截面密度、曲率尺寸以及端部重构。'], ...
             remeshAreaChange);

    end

end


fprintf('------------------------------------------------------------\n');



%% ========================================================================
% STEP 6 / 9
% 建立表面参数坐标
%
% 这一阶段建立：
%
%       x
%       r
%       theta
%       s
%
%
% 其中：
%
%       s = 从鼻部开始的几何子午弧长
%
%
% 注意：
%
% 这里的 s 还不是有攻角情况下的真实三维流线长度。
%% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 6 / 9 : Building body-fitted surface coordinates\n');
fprintf('============================================================\n');


mesh = ...
    buildSurfaceCoordinates( ...
        mesh, ...
        cfg);



%% ========================================================================
% 17. 检查 surfaceCoord 是否成功
% ========================================================================

if ~isfield(mesh,'surfaceCoord')

    error( ...
        'buildSurfaceCoordinates 没有生成 mesh.surfaceCoord。');

end


requiredSurfaceFields = { ...
    'nodeX', ...
    'nodeR', ...
    'nodeTheta', ...
    'nodeS', ...
    'faceX', ...
    'faceR', ...
    'faceTheta', ...
    'faceS'};


for i = 1:numel(requiredSurfaceFields)

    if ~isfield( ...
            mesh.surfaceCoord, ...
            requiredSurfaceFields{i})

        error( ...
            'mesh.surfaceCoord.%s 不存在。', ...
            requiredSurfaceFields{i});

    end

end



%% ========================================================================
% 18. 输出表面坐标统计
% ========================================================================

fprintf('\n');

fprintf('------------------------------------------------------------\n');
fprintf(' Surface Coordinate Summary\n');
fprintf('------------------------------------------------------------\n');


fprintf('Axial sections           : %d\n', ...
    mesh.surfaceCoord.numberSections);


fprintf('Common circum points     : %d\n', ...
    mesh.surfaceCoord.numberCommonCircumPoints);


fprintf('\n');


fprintf('Mean meridional length   : %.6f m\n', ...
    mesh.surfaceCoord.meanMeridionalLength);


fprintf('Minimum meridional length: %.6f m\n', ...
    mesh.surfaceCoord.minimumMeridionalLength);


fprintf('Maximum meridional length: %.6f m\n', ...
    mesh.surfaceCoord.maximumMeridionalLength);


fprintf('\n');


fprintf('Surface coordinate s     : Nose -> Aft\n');

fprintf('theta = 0                : +Y\n');

fprintf('theta = pi/2             : +Z\n');


fprintf('------------------------------------------------------------\n');



%% Local geometry on the accepted final mesh
mesh = analyzeLocalGeometry(mesh,cfg);


%% ========================================================================
% STEP 7 / 9
% Flight + Freestream GUI
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 7 / 9 : Flight and freestream definition\n');
fprintf('============================================================\n');


% 不适用压力模型的工况留在当前框架内修改，不留下失效页面。
while true
    setAeroHeatStage(appFigure,2);
[flight, ...
 freestream, ...
 flightAction] = ...
    createFlightCondition( ...
        cfg, ...
        mesh, appHost);
    if strcmpi(flightAction,'cancel') || freestream.Mach >= 1, break; end
    cfg.flight = flight;
    uialert(appFigure,'当前压力模型要求 Mach >= 1，请修改工况或取消。','模型适用范围');
end



%% ========================================================================
% 19. 用户取消飞行条件
% ========================================================================

if strcmpi( ...
        flightAction, ...
        'cancel')

    fprintf('\n');

    fprintf('============================================================\n');
    fprintf(' Flight-condition setup cancelled.\n');
    fprintf(' Main program terminated before aerodynamic calculation.\n');
    fprintf('============================================================\n');

    fprintf('\n');


    if isgraphics(appFigure), delete(appFigure); end
    return;

end



%% ========================================================================
% STEP 8 / 9
% 保存物理状态
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 8 / 9 : Storing atmospheric and freestream state\n');
fprintf('============================================================\n');


atmosphere = ...
    flight.atmosphere;



%% ========================================================================
% 20. 保存最终飞行输入到 cfg
% ========================================================================

cfg.flight.altitude = ...
    flight.altitude;


cfg.flight.Mach = ...
    flight.Mach;


cfg.flight.alpha = ...
    flight.alpha;


cfg.flight.beta = ...
    flight.beta;


cfg.flight.Twall = ...
    flight.Twall;



%% ========================================================================
% 21. 飞行 / 来流报告
% ========================================================================

fprintf('\n');

fprintf('------------------------------------------------------------\n');
fprintf(' Flight / Freestream Report\n');
fprintf('------------------------------------------------------------\n');


fprintf('Geometric altitude       : %.3f km\n', ...
    flight.altitude/1000);


fprintf('Geopotential altitude    : %.3f km\n', ...
    atmosphere.geopotentialAltitude/1000);


fprintf('\n');


fprintf('Mach number              : %.4f\n', ...
    flight.Mach);


fprintf('Angle of attack          : %.4f deg\n', ...
    flight.alpha);


fprintf('Sideslip angle           : %.4f deg\n', ...
    flight.beta);


fprintf('Wall temperature         : %.3f K\n', ...
    flight.Twall);


fprintf('\n');


fprintf('T infinity               : %.3f K\n', ...
    atmosphere.temperature);


fprintf('p infinity               : %.6e Pa\n', ...
    atmosphere.pressure);


fprintf('rho infinity             : %.6e kg/m^3\n', ...
    atmosphere.density);


fprintf('Speed of sound           : %.3f m/s\n', ...
    atmosphere.speedOfSound);


fprintf('Dynamic viscosity        : %.6e Pa*s\n', ...
    atmosphere.dynamicViscosity);


fprintf('\n');


fprintf('V infinity               : %.3f m/s\n', ...
    freestream.speed);


fprintf('q infinity               : %.6e Pa\n', ...
    freestream.dynamicPressure);


fprintf('Re / m                   : %.6e 1/m\n', ...
    freestream.reynoldsPerMeter);


fprintf('T0 infinity              : %.3f K\n', ...
    freestream.totalTemperature);


fprintf('\n');


fprintf('Vehicle direction        : [% .5f % .5f % .5f]\n', ...
    freestream.unitFlightDirection(1), ...
    freestream.unitFlightDirection(2), ...
    freestream.unitFlightDirection(3));


fprintf('Air flow direction       : [% .5f % .5f % .5f]\n', ...
    freestream.unitFlowDirection(1), ...
    freestream.unitFlowDirection(2), ...
    freestream.unitFlowDirection(3));


fprintf('------------------------------------------------------------\n');



%% 气动热区域分类：在用户确认飞行条件后，根据当前来流标记各面元
mesh = classifyAerothermalRegion(mesh,freestream,cfg);

%% 压力计算与可视化：飞行条件确认后计算，图窗关闭不影响工作区结果
if freestream.Mach >= 1
    mesh = modifiedNewtonianPressure(mesh,freestream,cfg);
    mesh = integrateAerodynamicLoads(mesh,cfg);
    mesh = findStagnationRegion(mesh,freestream,cfg);
    mesh = solveExternalFlow(mesh,freestream,cfg);
    mesh = marchAxisymmetricShockExpansion(mesh,freestream,cfg);
    mesh = traceSurfaceStreamlines(mesh,cfg);
    mesh = calcSurfaceHeating(mesh,freestream,cfg);
    fprintf('External-flow method     : %s\n',mesh.boundaryLayerEdge.method);
    if mesh.stagnationShock.available
        fprintf('Stagnation shock method  : %s\n', ...
            mesh.stagnationShock.method);
        fprintf('Post-shock p / T / Mach  : %.6e Pa / %.3f K / %.4f\n', ...
            mesh.stagnationShock.postShockPressure, ...
            mesh.stagnationShock.postShockTemperature, ...
            mesh.stagnationShock.postShockMach);
        fprintf('Shock total-p ratio      : %.6f\n', ...
            mesh.stagnationShock.totalPressureRatio);
    else
        fprintf('Stagnation shock status  : %s\n', ...
            mesh.stagnationShock.status);
    end
    if mesh.bowShockPreview.available
        fprintf('Bow-shock preview method: %s\n', ...
            mesh.bowShockPreview.method);
        fprintf('Estimated stand-off      : %.6e m (empirical)\n', ...
            mesh.bowShockPreview.standOffDistance);
    else
        fprintf('Bow-shock preview status: %s\n', ...
            mesh.bowShockPreview.status);
    end
    fprintf('External-flow valid faces: %d / %d\n', ...
        mesh.boundaryLayerEdge.summary.validFaces,size(mesh.F,1));
    fprintf('  tangent-cone faces     : %d\n',mesh.boundaryLayerEdge.summary.coneFaces);
    fprintf('  tangent-wedge faces    : %d\n',mesh.boundaryLayerEdge.summary.wedgeFaces);
    fprintf('  expansion faces        : %d\n',mesh.boundaryLayerEdge.summary.expansionFaces);
    fprintf('  detached-closure faces : %d (low-order)\n', ...
        mesh.boundaryLayerEdge.summary.detachedClosureFaces);
    fprintf('Axial shock-expansion    : %s\n',mesh.axisymmetricMarch.status);
    if mesh.axisymmetricMarch.applied
        fprintf('  marched faces          : %d\n', ...
            mesh.axisymmetricMarch.summary.marchedFaces);
    end
    angles=mesh.boundaryLayerEdge.shockAngleDeg;
    angles=angles(isfinite(angles));
    if ~isempty(angles)
        fprintf('Local shock angle range  : %.3f ~ %.3f deg\n', ...
            min(angles),max(angles));
    end
    fprintf('Streamline method        : %s\n',mesh.surfaceStreamlines.method);
    fprintf('Streamline valid faces   : %d / %d\n', ...
        mesh.surfaceStreamlines.summary.connectedFaces,size(mesh.F,1));
    fprintf('Boundary-layer method   : %s\n',mesh.boundaryLayer.method);
    fprintf('Laminar BL valid faces   : %d / %d\n', ...
        mesh.boundaryLayer.summary.validFaces,size(mesh.F,1));
    fprintf('Turbulent BL method     : %s\n',mesh.heating.turbulentMethod);
    if isstruct(mesh.heating.turbulentBoundaryLayer)
        turbulentCount=mesh.heating.turbulentBoundaryLayer.summary;
        fprintf('Turbulent BL sources    : %d handoff, %d integrated, %d local, %d fallback\n', ...
            turbulentCount.laminarHandoffFaces, ...
            turbulentCount.integratedFaces, ...
            turbulentCount.localStartupFaces, ...
            turbulentCount.flatPlateFallbackFaces);
        fprintf('Turbulent separation    : %d onset, %d including downstream; numerical gaps %d\n', ...
            turbulentCount.separationOnsetFaces,turbulentCount.separationUnresolvedFaces, ...
            turbulentCount.numericalUnresolvedFaces);
        fprintf('Turbulent domain blend  : %d faces (explicit engineering blend)\n', ...
            turbulentCount.blendedFaces);
    end
    fprintf('Edge-state method        : %s\n',mesh.boundaryLayerEdge.method);
    fprintf('Heating method           : %s\n',mesh.heating.method);
    fprintf('Air thermal properties   : %s\n',mesh.heating.thermalPropertyScope);
    fprintf('Continuum screen        : %s\n',mesh.continuumRegime.method);
    fprintf('Kn infinity / limit     : %.4g / %.4g (L=%.4g m)\n', ...
        mesh.continuumRegime.freestreamKnudsen, ...
        mesh.continuumRegime.maxContinuumKnudsen, ...
        mesh.continuumRegime.characteristicLength);
    fprintf('Continuum unresolved    : %d direct + %d downstream faces\n', ...
        mesh.continuumRegime.summary.directUnresolvedFaces, ...
        mesh.continuumRegime.summary.downstreamUnresolvedFaces);
    fprintf('Transition method        : %s\n',mesh.heating.transitionMethod);
    fprintf('Transition onset faces   : %d\n', ...
        mesh.heating.transition.summary.onsetFaces);
    fprintf('Inherited turbulent faces: %d\n', ...
        mesh.heating.transition.summary.inheritedTurbulentFaces);
    fprintf('Transition through heat gaps: %d carrier / %d downstream faces\n', ...
        mesh.heating.transition.summary.turbulentCarrierFaces, ...
        mesh.heating.transition.summary.inheritedThroughUnavailableFaces);
    fprintf('Calculated heat faces    : %d / %d\n', ...
        mesh.heating.summary.calculatedFaces,size(mesh.F,1));
    fprintf('Calculated heat area     : %.2f %% of mesh surface\n', ...
        100*mesh.heatingCoverage.summary.computedAreaFraction);
    fprintf('Nose extension method    : %s\n',mesh.heating.noseMethod);
    fprintf('Nose extension faces     : %d (low-order estimate)\n', ...
        mesh.heating.summary.noseFaces);
    fprintf('Detached-closure heat faces: %d (low-order estimate)\n', ...
        mesh.heating.summary.detachedClosureHeatFaces);
    fprintf('Stagnation method        : %s\n',mesh.heating.stagnationMethod);
    fprintf('Stagnation status        : %s\n',mesh.heating.stagnationStatus);
    fprintf('Nose radius source       : %s\n',mesh.stagnationHeating.radiusSource);
    if mesh.stagnationHeating.available
        fprintf('Nose radius              : %.6e m\n', ...
            mesh.stagnationHeating.radius);
        fprintf('Stagnation heat flux     : %.6e W/m^2\n', ...
            mesh.stagnationHeating.heatFlux);
        fprintf('Stagnation incidence     : %.3f deg\n', ...
            mesh.stagnationHeating.incidenceAngleDeg);
        fprintf('Stagnation position      : [% .6e % .6e % .6e] m\n', ...
            mesh.stagnationHeating.position);
    end
    setAeroHeatStage(appFigure,3);
    pressureFigure = pressureGUI(mesh,cfg,appHost);
else
    warning('RocketAeroHeat:Pressure:Subsonic', ...
        '当前为亚声速，已保留前处理结果；修正牛顿压力模块要求 Mach >= 1。');
end

%% ========================================================================
% STEP 9 / 9
% 当前前处理阶段完成
% ========================================================================

fprintf('\n');

fprintf('============================================================\n');
fprintf('Step 9 / 9 : Preprocessing completed\n');
fprintf('============================================================\n');


fprintf('\n');


fprintf('Important variables:\n');

fprintf('\n');


fprintf('    rawMeshOriginal\n');
fprintf('        Original STL coordinate system\n');

fprintf('\n');


fprintf('    rawMesh\n');
fprintf('        Geometry aligned to body X axis\n');

fprintf('\n');


fprintf('    mesh\n');
fprintf('        Final computational surface mesh\n');

fprintf('\n');


fprintf('    mesh.param\n');
fprintf('        Native axial-section mesh topology\n');

fprintf('\n');


fprintf('    mesh.surfaceCoord\n');
fprintf('        x / r / theta / s body-fitted coordinates\n');

fprintf('\n');


fprintf('    quality\n');
fprintf('        Final triangle quality\n');

fprintf('\n');


fprintf('    cfg\n');
fprintf('        Full program configuration\n');

fprintf('\n');


fprintf('    flight\n');
fprintf('        User-defined flight condition\n');

fprintf('\n');


fprintf('    atmosphere\n');
fprintf('        Standard atmospheric state\n');

fprintf('\n');


fprintf('    freestream\n');
fprintf('        3-D freestream state and direction\n');



%% ========================================================================
% 22. 后续开发提示
% ========================================================================

fprintf('\n');

fprintf('Coordinate convention used by later modules:\n');

fprintf('\n');

fprintf('    Body axis             = X\n');

fprintf('    Nose                  = +X\n');

fprintf('    +Y                    = lateral\n');

fprintf('    +Z                    = upward\n');


fprintf('\n');


fprintf('Surface coordinates available:\n');

fprintf('\n');

fprintf('    mesh.surfaceCoord.faceX\n');

fprintf('    mesh.surfaceCoord.faceR\n');

fprintf('    mesh.surfaceCoord.faceTheta\n');

fprintf('    mesh.surfaceCoord.faceS\n');


fprintf('\n');


fprintf('Current engineering modules:\n');

fprintf('\n');

fprintf('    analyzeLocalGeometry.m : completed (mesh.localGeometry)\n');

fprintf('    classifyAerothermalRegion.m : completed (mesh.aerothermalRegion)\n');


fprintf('\n');


fprintf('Aerodynamic models:\n');

fprintf('\n');

fprintf('    calcCpMax.m / modifiedNewtonianPressure.m : integrated\n');

fprintf('    pressureGUI.m : pressure and heating viewer\n');

fprintf('\n');


fprintf('Aerothermal models:\n');

fprintf('\n');

fprintf('    findStagnationRegion.m : candidate regions available in mesh.stagnation\n');

fprintf('    calcBoundaryLayerEdgeState.m : engineering edge-state estimate\n');

fprintf('    calcBoundaryLayerState.m : streamline Thwaites-type laminar integral\n');
fprintf('    calcTurbulentBoundaryLayerState.m : streamline turbulent momentum integral\n');
fprintf('    calcSurfaceHeating.m : integrated laminar/turbulent heat flux\n');
fprintf('    assessContinuumRegime.m : Knudsen applicability screening\n');
fprintf('    calcBluntNoseHeating.m : offset-stagnation spherical nose extension\n');

fprintf('    estimateNoseRadius.m / calcStagnationHeating.m : validated nose radius + Sutton-Graves\n');

fprintf('    solveObliqueShock.m / solveConicalShock.m : attached-shock solvers\n');
    fprintf('    solveNormalShock.m : local blunt-nose stagnation shock state\n');
    fprintf('    solveThermallyPerfectNormalShock.m : NASA7 local nose shock\n');
    fprintf('    solveThermallyPerfectConicalShock.m : NASA7 zero-alpha cone shock\n');
    fprintf('    solveThermallyPerfectObliqueShock.m : NASA7 local wedge shock\n');
    fprintf('    solveThermallyPerfectPrandtlMeyerExpansion.m : NASA7 local expansion\n');
    fprintf('    recoverThermallyPerfectEdgeState.m : NASA7 isentropic edge state\n');
    fprintf('    airThermalProperties.m : dry-air cp / h / entropy\n');
fprintf('    estimateBowShockSurface.m : Billig blunt-nose empirical shock preview\n');

fprintf('    solveExternalFlow.m : built-in engineering external flow\n');
fprintf('    marchAxisymmetricShockExpansion.m : zero-alpha axial wave march\n');

fprintf('    importExternalFlowState.m : validated CFD edge-state import\n');

fprintf('    traceSurfaceStreamlines.m : velocity-directed surface distances\n');
fprintf('    traceSurfacePath.m : selected three-dimensional discrete path\n');
fprintf('    calcFlightHeatingHistory.m : built-in external flow at each flight station\n');
fprintf('    refineFlightStations.m : linear input interpolation and bounded time steps\n');
fprintf('    calcWallThermalHistory.m : variable-condition wall temperature history\n');
fprintf('    calcWallThermalResponse.m : fixed-condition wall temperature history\n');
fprintf('    calcAerothermalExposure.m : integrated heat dose and energy check\n');
fprintf('    analyzeHeatingCoverage.m : heat-result area and unresolved-reason audit\n');
fprintf('    pressureGUI.m : pressure, heating and thermal-history pages\n');


fprintf('\n');

fprintf('============================================================\n');
fprintf(' Rocket AeroHeat preprocessing completed successfully.\n');
fprintf('============================================================\n');

fprintf('\n');



%% ========================================================================
% LOCAL FUNCTION
%
% 检查 Boundary Edge / Non-manifold Edge
%% ========================================================================

function [boundaryEdges,nonManifoldEdges] = ...
    localCheckTopology(F)


%% ------------------------------------------------------------------------
% 所有三角形边
%% ------------------------------------------------------------------------

edges = [
    F(:,[1 2]);
    F(:,[2 3]);
    F(:,[3 1])
];


%% ------------------------------------------------------------------------
% 转换成无向边
%% ------------------------------------------------------------------------

edges = ...
    sort( ...
        edges, ...
        2);


%% ------------------------------------------------------------------------
% 公共边分组
%% ------------------------------------------------------------------------

[~,~,group] = ...
    unique( ...
        edges, ...
        'rows');


%% ------------------------------------------------------------------------
% 每条边出现次数
%% ------------------------------------------------------------------------

edgeCount = ...
    accumarray( ...
        group, ...
        1);


%% ------------------------------------------------------------------------
% 出现一次：
%
% Boundary Edge
%% ------------------------------------------------------------------------

boundaryEdges = ...
    sum( ...
        edgeCount == 1);


%% ------------------------------------------------------------------------
% 出现 >2 次：
%
% Non-manifold Edge
%% ------------------------------------------------------------------------

nonManifoldEdges = ...
    sum( ...
        edgeCount > 2);


end
