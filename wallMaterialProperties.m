function p = wallMaterialProperties(material,temperature)
%WALLMATERIALPROPERTIES 常物性或用户温度表的材料比热/导热系数。
% 模块一：校验 propertyTable.temperature/specificHeat/thermalConductivity。
% 模块二：在表内分段线性插值 cp [J/(kg K)] 和 k [W/(m K)]。
% 表外禁止外推；NaN 温度表示未解面元，返回 NaN，不补成初始温度。
% 密度保持常数；不含相变潜热、热膨胀和材料化学反应。
validateattributes(temperature,{'numeric'},{'real'});
finite=isfinite(temperature);
if any(isinf(temperature),'all') || any(temperature(finite)<=0)
    error('RocketAeroHeat:Material:Temperature','材料温度必须为正有限值或 NaN。');
end
if ~isfield(material,'propertyTable') || isempty(material.propertyTable)
    validateattributes(material.specificHeat,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(material.thermalConductivity,{'numeric'},{'real','finite','scalar','positive'});
    cp=material.specificHeat+zeros(size(temperature));
    k=material.thermalConductivity+zeros(size(temperature));
    cp(~finite)=NaN; k(~finite)=NaN;
    p=struct('specificHeat',cp,'thermalConductivity',k, ...
        'temperatureRange',[0 Inf],'method','用户给定常比热/常导热系数');
    return;
end
tab=material.propertyTable;
keys={'temperature','specificHeat','thermalConductivity'};
for i=1:numel(keys)
    if ~isfield(tab,keys{i})
        error('RocketAeroHeat:Material:Table','物性表缺少 %s。',keys{i});
    end
    validateattributes(tab.(keys{i}),{'numeric'}, ...
        {'real','finite','vector','nonempty','positive'});
end
T=double(tab.temperature(:)); cp=double(tab.specificHeat(:)); k=double(tab.thermalConductivity(:));
if numel(T)<2 || any(diff(T)<=0) || numel(cp)~=numel(T) || numel(k)~=numel(T)
    error('RocketAeroHeat:Material:Table','物性表至少两行，温度严格递增，三列等长。');
end
if any(temperature(finite)<T(1) | temperature(finite)>T(end))
    error('RocketAeroHeat:Material:Range', ...
        '材料温度超出用户物性表 %.6g~%.6g K；禁止外推，请补充有效数据或缩短计算时间。',T(1),T(end));
end
p=struct('specificHeat',reshape(interp1(T,cp,temperature(:),'linear'),size(temperature)), ...
    'thermalConductivity',reshape(interp1(T,k,temperature(:),'linear'),size(temperature)), ...
    'temperatureRange',[T(1) T(end)], ...
    'method','用户 cp(T)、k(T) 表的分段线性插值；表外不外推');
end
