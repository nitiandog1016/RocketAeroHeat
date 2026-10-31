function [options,mode]=su2WorkbenchMeshOptions(cfg)
%SU2WORKBENCHMESHOPTIONS 将可见网格模式转换为已测试的受限原生箭体参数。
% 模块一：新配置显式选择法向偏置层；旧快照无模式字段时保持原射线网格。
% 模块二：只提供经过几何/导入检查的参数组合，不赋予近壁或物理精度资格。
mode='starRay';
if isfield(cfg,'su2') && isfield(cfg.su2,'volumeMode'),mode=char(cfg.su2.volumeMode);end
options=struct('radialLayers',8,'maxWallAngleStepDeg',10,'baseRadialSegments',4);
switch mode
    case 'normalLayers'
        options.radialLayers=40;options.firstRadialSpacing=2e-6;
        options.nearWall=struct('normalLayers',24,'normalGrowth',1.25);
    case 'starRay'
        % 旧模式参数保持不变，用于历史快照和网格对照。
    case 'hybridTriangles'
        % 原三角形近壁层+Gmsh外域，不要求轴对称母线；仍须逐模型几何审核。
        options=struct;
        if isfield(cfg.su2,'hybridMesh'),options=cfg.su2.hybridMesh;end
        if isfield(cfg.su2,'pythonExecutable'),options.pythonExecutable=cfg.su2.pythonExecutable;end
    otherwise
        error('RocketAeroHeat:SU2:VolumeMode','不支持的体网格模式：%s',mode);
end
end
