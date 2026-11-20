function [b,fs]=thermistorWorkbenchBoundary(b)
%THERMISTORWORKBENCHBOUNDARY 给出两档周向网共同的等温试件工况，拒绝错配边界。
reference=struct('speed',3,'airPressurePa',101325,'airTemperatureK',293.15, ...
    'wallTemperatureK',305.15,'alphaDeg',0,'betaDeg',0);
if nargin==0,b=reference;end
assert(isequal(b,reference),'RAH:Bead:Boundary','此热敏珠算例仅支持登记的 3 m/s 等温边界。');
fs=struct('speed',b.speed,'pressure',b.airPressurePa,'temperature',b.airTemperatureK, ...
    'alpha',0,'beta',0,'unitFlowDirection',[-1 0 0],'velocityMps',[-b.speed 0 0]);
end
