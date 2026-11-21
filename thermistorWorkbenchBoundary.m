function [b,fs]=thermistorWorkbenchBoundary(b,angularCells,cfg)
%THERMISTORWORKBENCHBOUNDARY 核对空气来流、互斥热边界及所选网格的净功率上限。
reference=struct('speed',3,'airPressurePa',101325,'airTemperatureK',293.15, ...
    'wallTemperatureK',305.15,'alphaDeg',0,'betaDeg',0);
if nargin<2,angularCells=2;end
if nargin<3,cfg=struct();end
if angularCells==48
    % 全周启动工况的壁温和来流与构造初场一致，限定终点由两个物理半步确定。
    reference.speed=6;reference.wallTemperatureK=303.15;
    assert(isequal(b,reference),'RAH:Bead:Boundary','56 层全周网使用 6 m/s、293.15 K 空气和 303.15 K 给定壁温。');
end
if nargin==0,b=reference;end
assert(isstruct(b) && isscalar(b) && xor(isfield(b,'wallTemperatureK'),isfield(b,'netHeatingPowerW')), ...
    'RAH:Bead:Boundary','请选择一项热边界：均匀壁温或整珠净加热功率。');
key='wallTemperatureK';powerMode=isfield(b,'netHeatingPowerW');
if powerMode,key='netHeatingPowerW';end
if isfield(cfg,'thermistorVariant') && any(strcmp(cfg.thermistorVariant,{'leads','midbody'}))
    % 第二、第三档带线网支持两种热边界；其余混合网须使用原生给定壁温。
    adjustable=strcmp(cfg.thermistorVariant,'leads') && isfield(cfg,'thermistorRefineLeads') && ...
        ismember(cfg.thermistorRefineLeads,[1 2]);
    assert(angularCells==4 && (adjustable || isequal(b,reference)), ...
        'RAH:Lead:Boundary','请选择与当前试件体网匹配的热边界。');
end
assert(isequal(rmfield(b,key),rmfield(reference,'wallTemperatureK')), ...
    'RAH:Bead:Boundary','此算例使用登记的 3 m/s 空气条件。');
value=b.(key);
assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value), ...
    'RAH:Bead:Boundary','热边界输入须为有限实数标量。');
value=double(value);
if powerMode
    limit=thermistorNetPowerLimit(angularCells,cfg);
    assert(value>=0 && value<=limit,'RAH:Bead:Boundary', ...
        '该网格的净加热功率须在 0 至 %.15g W 之间。',limit);
else
    assert(value>=reference.airTemperatureK && value<=reference.wallTemperatureK, ...
        'RAH:Bead:Boundary','给定壁温须为 293.15 至 305.15 K 的有限实数。');
end
b.(key)=value;
fs=struct('speed',b.speed,'pressure',b.airPressurePa,'temperature',b.airTemperatureK, ...
    'alpha',0,'beta',0,'unitFlowDirection',[-1 0 0],'velocityMps',[-b.speed 0 0]);
end
