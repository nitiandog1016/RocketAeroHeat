function volume=buildSU2HybridVolume(surface,folder,options)
%BUILDSU2HYBRIDVOLUME 由完整原表面自动构造棱柱近壁层及四面体外域。
% 模块一：原V/F和显式网格设置写入小型输入，不采用外部CFD场。
% 模块二：调用固定Gmsh子进程，限制2线程/720秒；失败证据保留。
% 模块三：独立读取并审核全部体单元、界面与原面映射，再交给SU2。
if nargin<3,options=struct;end
defaults=struct('layerHeights',[0 cumsum(2e-5*1.4.^(0:5))], ...
    'downstreamLengths',2,'upstreamLengths',1,'lateralLengths',1, ...
    'outerMinSize',.02,'farSize',.35,'pythonExecutable','python','maximumSeconds',120);
assert(all(ismember(fieldnames(options),fieldnames(defaults))),'不支持的混合网格参数。');
for name=fieldnames(defaults)',k=name{1};if ~isfield(options,k),options.(k)=defaults.(k);end,end
for name={'downstreamLengths','upstreamLengths','lateralLengths','outerMinSize','farSize','maximumSeconds'}
    validateattributes(options.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
assert(options.maximumSeconds<=720 && options.farSize>=options.outerMinSize);
folder=char(java.io.File(folder).getCanonicalPath());assert(~isfolder(folder),'不能覆盖已有网格任务。');
parent=fileparts(folder);if ~isfolder(parent),mkdir(parent);end
input=[folder '_input.mat'];assert(~isfile(input),'不能覆盖已有输入。');
V=surface.V;F=surface.F;python=options.pythonExecutable;seconds=options.maximumSeconds;
options=rmfield(options,{'pythonExecutable','maximumSeconds'});optionsJSON=jsonencode(options);
save(input,'V','F','optionsJSON','-v7');
script=fullfile(fileparts(mfilename('fullpath')),'scripts','meshSU2HybridVolume.py');
for value={python,script,input,folder}
    assert(~any(ismember(value{1},char([34 10 13 37]))),'执行路径含不支持字符。');
end
[code,message]=system(sprintf('"%s" "%s" --input "%s" --folder "%s" --maximum-seconds %.9g',python,script,input,folder,seconds));
assert(code==0,'RocketAeroHeat:SU2:HybridMesher','体网格构造失败，原证据保留在%s；%s',folder,message);
volume=readSU2HybridVolume(folder,surface);
end
