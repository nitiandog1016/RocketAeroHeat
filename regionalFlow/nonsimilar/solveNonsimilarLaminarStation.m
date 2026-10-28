function result=solveNonsimilarLaminarStation(initial,history,gas,domain,options)
%SOLVENONSIMILARLAMINARSTATION 有明确适用声明的非相似层流单站研究接口。
% 模块一：核对调用者声明，限制在连续流、轴对称、零攻角/侧滑、层流、
% 无逆流/分离及已检查薄层的研究条件；声明并不等于本函数自动验证物理。
% 模块二：调用既有物理n隐式推进，不改离散方程、物性、阈值或上游剖面。
% 模块三：保留入口/完整壁温历史身份与数值状态；固定physicalSolution和
% eligibleForSurfaceCoverage为false，数值收敛不能自动恢复生产解面积。
% 使用：addpath(projectRoot,fullfile(projectRoot,'regionalFlow','nonsimilar'))。
% domain必填字段见本文件required；入口必须是同物理n上的真实u/T剖面，
% 不支持驻点s=0、三维横流、分离再附、湍流、稀薄流或制造源。
% 物理适用性、逐站守恒账及公开数据验证仍需调用方独立完成。
if nargin<5,options=struct;end
if ~isstruct(domain)||~isscalar(domain)
    error('RocketAeroHeat:NonSimilar:Domain','domain必须是单一结构体适用声明。');
end
required={'axisymmetric','angleOfAttackDeg','sideslipDeg','flowRegime', ...
    'layerState','crossflowPresent','separationPresent','thinLayerVerified', ...
    'inletIdentity','wallTemperatureHistoryIdentity'};
for j=1:numel(required)
    if ~isfield(domain,required{j})
        error('RocketAeroHeat:NonSimilar:Domain','缺少domain.%s声明。',required{j});
    end
end
unknown=setdiff(fieldnames(domain),required);
if ~isempty(unknown)
    error('RocketAeroHeat:NonSimilar:Domain','未知domain字段%s；避免静默忽略边界要求。',unknown{1});
end
for key={'axisymmetric','crossflowPresent','separationPresent','thinLayerVerified'}
    z=domain.(key{1});
    if ~islogical(z)||~isscalar(z)
        error('RocketAeroHeat:NonSimilar:Domain','domain.%s必须为逻辑标量。',key{1});
    end
end
for key={'angleOfAttackDeg','sideslipDeg'}
    z=domain.(key{1});
    if ~isnumeric(z)||~isreal(z)||~isscalar(z)||~isfinite(z)
        error('RocketAeroHeat:NonSimilar:Domain','角度必须是有限实数标量，单位deg。');
    end
end
for key={'flowRegime','layerState','inletIdentity','wallTemperatureHistoryIdentity'}
    z=domain.(key{1});
    if ~(ischar(z)&&isrow(z) || isstring(z)&&isscalar(z)) || ...
            ismissing(string(z)) || strlength(strtrim(string(z)))==0
        error('RocketAeroHeat:NonSimilar:Domain','domain.%s必须为非空单一文本。',key{1});
    end
end
if ~domain.axisymmetric || domain.angleOfAttackDeg~=0 || domain.sideslipDeg~=0 || ...
        ~strcmp(domain.flowRegime,'continuum') || ~strcmp(domain.layerState,'laminar') || ...
        domain.crossflowPresent || domain.separationPresent || ~domain.thinLayerVerified
    error('RocketAeroHeat:NonSimilar:Domain', ...
        '本研究接口仅接受轴对称、严格零攻角/侧滑、连续流、层流、无横流/分离和已检查薄层的声明。');
end
% 模块二：不提供制造源参数；制造解只允许直接调用低层数值测试接口。
result=advanceNonsimilarLaminarStation(initial,history,gas,options);
% 模块三：区分调用者声明、离散收敛与物理区域有效性。
result.numericalConverged=result.converged;
result.physicalSolution=false;
result.eligibleForSurfaceCoverage=false;
result.domainDeclaration=domain;
result.domainDeclarationAutomaticallyVerified=false;
result.domainDeclarationScope='调用者适用声明；本接口没有自动验证轴对称/薄层/层流/无分离等物理条件';
result.moduleIdentity='v61受限条件非相似层流研究模块';
result.scope=['条件研究结果；数值收敛仅指当前离散方程。本接口不证明物理适用、' ...
    '驻点启动、精确能量守恒、全程稳定性、实验精度或新增生产面积。'];
end
