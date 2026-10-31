function folder=setupCutCellFlow
%SETUPCUTFLOW 加入完整曲面控制体实验库及共享特征边界模块。
% 模块用途：在任意当前工作目录下从本文件位置解析依赖，不依赖outputs。
% 只设置路径，不自动计算、不接管生产工程方法，也不更改结果覆盖身份。
parent=fileparts(mfilename('fullpath'));folder=fullfile(parent,'cutCell');
assert(isfolder(folder),'CutLibrary:Missing','缺失完整曲面数值库目录。');
addpath(parent,'-begin');addpath(folder,'-begin');
end
