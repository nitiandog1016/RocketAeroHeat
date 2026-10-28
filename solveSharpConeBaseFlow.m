function result = solveSharpConeBaseFlow(fs,geometry,Twall,cfg)
%SOLVESHARPCONEBASEFLOW 内部计算直圆锥肩部状态，再计算底部工程关联。
% 模块用途：不导入 CFD；使用本项目 Taylor-Maccoll 解生成直圆锥
% 肩部 p/T/rho/U/M，以锥母线长度生成 Re_e,s，然后调用有来源的
% solveAxisymmetricBaseCorrelation。仅用于几何已明确的尖头直圆锥，
% 不将其外流替代球鼻、锥柱、尾锥或非零攻角的真实肩部状态。
% geometry 需 baseRadius, halfAngleDeg, boundaryLayerRegime,
% hasRearSupport, hasGasInjection, source。底热为代表量，没有全底面
% 径向热分布；给定壁温、总焓参考和未解外缘静温分开存储。

%% 1. 声明尖头直锥几何；拒绝冒用钝锥肩部状态
if nargin<4, cfg=struct; end
if ~isfield(geometry,'halfAngleDeg') || ~isfield(geometry,'baseRadius')
    error('RocketAeroHeat:SharpConeBase:Geometry','需要半锥角与底半径。');
end
if isfield(geometry,'noseBluntnessRatio') && geometry.noseBluntnessRatio~=0
    error('RocketAeroHeat:SharpConeBase:Bluntness', ...
        '本内部外流封装只适用尖头直锥；钝锥需要其自己的肩部解。');
end
if isfield(geometry,'family') && ~strcmp(geometry.family,'cone')
    error('RocketAeroHeat:SharpConeBase:Family','本封装只适用直圆锥。');
end
geometry.family='cone'; geometry.noseBluntnessRatio=0;
angle=geometry.halfAngleDeg;
validateattributes(angle,{'numeric'},{'real','finite','scalar','>',0,'<',90});
validateattributes(geometry.baseRadius,{'numeric'}, ...
    {'real','finite','scalar','positive'});
% 下层关联在物理门限不满足时返回 NaN；此处仍需要可检查的状态。
state=solveConicalShock(fs,angle);
if ~state.attached
    result=struct('available',false,'pressureAvailable',false, ...
        'heatFluxAvailable',false,'basePressure',NaN, ...
        'heatFluxRepresentative',NaN,'eligibleForSurfaceCoverage',false, ...
        'status',['内部直圆锥肩部无附体弱解：' state.status], ...
        'internalShoulderState',state);
    return;
end

%% 2. 肩部外流和鼻端到肩部母线长度的逐输入来源
shoulder=struct;
for key={'pressure','temperature','density','speed','Mach'}
    shoulder.(key{1})=state.(key{1});
    shoulder.sources.(key{1})=['内部求解：' state.method '；尖头直圆锥锥面状态'];
end
shoulder.dynamicViscosity=1.458e-6*state.temperature^1.5/(state.temperature+110.4);
shoulder.distance=geometry.baseRadius/sind(angle);
shoulder.sources.dynamicViscosity='Sutherland 空气黏度，与生产输运实现一致';
shoulder.sources.distance='尖头直圆锥几何母线长度 rb/sin(theta_c)';

%% 3. 区域关联与代表底热；不把标量映射为未验证的底面热场
result=solveAxisymmetricBaseCorrelation(fs,shoulder,geometry,Twall,cfg);
result.internalShoulderState=state;
result.internalExternalFlowComputed=true;
end
