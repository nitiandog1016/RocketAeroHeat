function setupEngineeringShockLayer()
%SETUPENGINEERINGSHOCKLAYER 加载内部三维工程激波层模块及同物性接口。
% 模块用途：仅设置路径；不替换旧外流、网格、热流或正式资格状态。
root=fileparts(mfilename('fullpath'));
setupSurfacePressureFlow;
addpath(fullfile(root,'regionalFlow','shockLayer'));
end
