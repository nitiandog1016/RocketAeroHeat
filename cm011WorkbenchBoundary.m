function [b,fs]=cm011WorkbenchBoundary(b)
%CM011WORKBENCHBOUNDARY 返回已验证低速试件的边界，并核对界面实际选择的值。
reference=struct('speed',0.1,'alphaDeg',10,'betaDeg',0,'airPressurePa',101325, ...
    'airTemperatureK',300,'wallTemperatureK',300.1);
if nargin==0,b=reference;end
keys=fieldnames(reference);
assert(isstruct(b)&&isscalar(b)&&all(isfield(b,keys)),'RAH:CM011:Boundary','试件边界不完整。');
for k=1:numel(keys)
    key=keys{k};
    assert(isequal(b.(key),reference.(key)),'RAH:CM011:Boundary', ...
        '当前 CM-011 适配器支持 %s = %.9g，请核对边界。',key,reference.(key));
end
direction=[-cosd(b.alphaDeg) 0 sind(b.alphaDeg)];
fs=struct('speed',b.speed,'pressure',b.airPressurePa,'temperature',b.airTemperatureK, ...
    'alpha',b.alphaDeg,'beta',b.betaDeg,'unitFlowDirection',direction,'velocityMps',b.speed*direction);
end
