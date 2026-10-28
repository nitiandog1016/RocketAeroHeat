function [value,report,evidence]=solveRocketTemperatureActiveDirection(packed,factor,rhs,metric,rows,options)
%SOLVEROCKETTEMPERATUREACTIVEDIRECTION 从自由方向开始，逐次加入朝原温域外变化的界面行。
% 模块一：每个当前状态从零自由方向重算，上一状态的约束与导数不继承。
% 模块二：投影后若其他带内行转向界外，再加入并重算；当前状态内不删已激活行。
% 模块三：只是可行方向原型，不声称不等式最优解或非线性温域保证，原NS仍独立验收。
if nargin<6,options=struct;end
% 模块四：显式启用时，用本次刚计算的自由方向作前瞻选面；缺省路径和输出不变。
lookahead=[];
if isstruct(options)&&isscalar(options)&&isfield(options,'lookaheadPreview')
    lookahead=options.lookaheadPreview;options=rmfield(options,'lookaheadPreview');
    required={'Q','info','cut','gas','scale','targetAlpha'};allowed=[required {'selectionPolicy'}];
    assert(isstruct(lookahead)&&isscalar(lookahead)&&all(isfield(lookahead,required))&& ...
        all(ismember(fieldnames(lookahead),allowed)), ...
        'TemperatureLookahead:Input','前瞻请求必须明确同一当前状态、采样、几何、气体、缩放和目标步长。');
    assert(isequal(size(lookahead.Q),[packed.cellCount 5])&&rows.options.temperatureBandK==10, ...
        'TemperatureLookahead:Input','前瞻仅接受当前全部CV和原10K基础行，不继承其他筛选。');
    validateattributes(lookahead.targetAlpha,{'numeric'},{'scalar','real','finite','positive','<=',1});
end
defaults=struct('maximumIterations',128,'maximumActiveSolves',3,'checkpointFolder','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'TemperatureActive:Options','未知活跃方向选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(options.maximumActiveSolves,{'numeric'},{'scalar','integer','positive','<=',3});
assert(rows.cellCount==packed.cellCount&&all(ismember(rows.temperatureBoundSigns,[-1 1]))&& ...
    numel(rows.temperatureBoundSigns)==size(rows.H,1),'TemperatureActive:Rows','当前温度行或上/下界符号不符。');
if ~isempty(options.checkpointFolder),assert(isfolder(options.checkpointFolder));end
[free,freeDetail,freeCoordinates]=solveRocketAugmentedQRDirection(packed,factor,rhs,metric,struct('maximumIterations',options.maximumIterations));
lookaheadReport=[];
if ~isempty(lookahead)
    originalRows=rows;
    if ~isempty(options.checkpointFolder)
        % 前瞻身份检查若拒绝，仍保留同一当前系统的自由前缀，避免丢失昂贵计算。
        prefix=fullfile(options.checkpointFolder,'lookahead_free_prefix.mat');assert(~isfile(prefix));
        save(prefix,'free','freeDetail','freeCoordinates','originalRows','-v7');
    end
    selectionOptions=struct;
    if isfield(lookahead,'selectionPolicy'),selectionOptions.selectionPolicy=lookahead.selectionPolicy;end
    [rows,lookaheadReport]=selectRocketTemperatureLookaheadRows(lookahead.Q,lookahead.info,lookahead.cut, ...
        lookahead.gas,lookahead.scale,free,lookahead.targetAlpha,selectionOptions);
    [retained,locations]=ismember(originalRows.faceIDs,rows.faceIDs);
    assert(all(retained)&&isequal(originalRows.H,rows.H(locations,:))&& ...
        isequal(originalRows.temperatureBoundSigns,rows.temperatureBoundSigns(locations)), ...
        'TemperatureLookahead:Retention','调用者原10K条件未逐项保留，拒绝当前请求。');
end
[value,report,evidence]=finishRocketTemperatureActiveDirection(packed,factor,rhs,metric,rows,options,free,freeDetail,freeCoordinates);
if ~isempty(lookaheadReport),report.lookaheadSelection=lookaheadReport;end
end
