function base=readRocketNSProbeBase(file,stateLocation)
%READROCKETNSPROBEBASE 按显式root/result位置读取真实基态，不改写检查点历史。
% 模块一：默认root保留已有顶层Q/rate合同；result须明确指定最终接受状态。
% 模块二：只提取真实Q/rate/尺度/物性/边界，不制造新的求解器result结构。
assert(ischar(file)&&isfile(file)&&ischar(stateLocation)&&ismember(stateLocation,{'root','result'}), ...
    'SteadyGuess:BaseSource','基态来源必须为实际MAT文件与明确root/result位置。');
available=whos('-file',file);names={available.name};
assert(all(ismember({'gas','boundary'},names)),'SteadyGuess:BaseSource','来源缺少物性或边界。');
if strcmp(stateLocation,'result')
    assert(ismember('result',names),'SteadyGuess:BaseSource','来源没有实际result，不能从别处补造。');
    s=load(file,'result','gas','boundary');r=s.result;
    assert(all(isfield(r,{'Q','rate','residualScale'})),'SteadyGuess:BaseSource','实际result缺少完整末态/残差/尺度。');
    base=struct('Q',r.Q,'rate',r.rate,'scale',r.residualScale,'gas',s.gas,'boundary',s.boundary);
else
    assert(all(ismember({'Q','rate'},names)),'SteadyGuess:BaseSource','root来源必须实际保存Q和rate，不自动切换到result。');
    keep=intersect({'Q','rate','scale','result','gas','boundary'},names,'stable');s=load(file,keep{:});
    if isfield(s,'scale'),scale=s.scale;
    else,assert(isfield(s,'result')&&isfield(s.result,'residualScale'),'SteadyGuess:BaseSource','root来源没有明确尺度。');scale=s.result.residualScale;end
    base=struct('Q',s.Q,'rate',s.rate,'scale',scale,'gas',s.gas,'boundary',s.boundary);
end
assert(isequal(size(base.Q),size(base.rate))&&size(base.Q,2)==5&&numel(base.scale)==5&&all(isfinite(base.scale)&base.scale>0), ...
    'SteadyGuess:BaseSource','真实基态/残差/尺度尺寸非法。');
base.scale=reshape(base.scale,1,5);
end
