function setupShockFittedFlow()
%SETUPSHOCKFITTEDFLOW 加载局部三维激波拟合/ALE研究内核及同物性接口。
% 模块用途：只设置路径；不替换生产外流、网格或区域资格。
root=fileparts(mfilename('fullpath'));setupEngineeringShockLayer;
addpath(fullfile(root,'internalFlow'));setupCutCellFlow;
addpath(fullfile(root,'regionalFlow','shockFitting'));
end
