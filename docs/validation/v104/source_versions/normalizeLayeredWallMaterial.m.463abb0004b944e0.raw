function material=normalizeLayeredWallMaterial(material,initial)
%NORMALIZELAYEREDWALLMATERIAL 验证分层输入并生成兼容旧结构的参考量。
% 以下等效密度/比热/导热系数仅供输入记录和旧字段兼容；多层求解
% 与储热核对使用逐层物性，绝不按这些平均量求温度。
validateattributes(initial,{'numeric'},{'real','finite','scalar','positive'});
grid=buildLayeredWallGrid(material);
mass=0; capacity=0; thermalResistance=sum(grid.contactResistance);
for j=1:numel(grid.layers)
    layer=grid.layers{j}; p=wallMaterialProperties(layer,initial);
    localMass=layer.density*layer.thickness;
    mass=mass+localMass; capacity=capacity+localMass*p.specificHeat;
    thermalResistance=thermalResistance+layer.thickness/p.thermalConductivity;
end
material.layers=grid.layers;
material.contactResistance=grid.contactResistance;
material.thickness=grid.depth(end);
material.density=mass/material.thickness;
material.specificHeat=capacity/mass;
material.thermalConductivity=material.thickness/thermalResistance;
material.thermalNodes=grid.nodeCount;
% 顶层单材料表与逐层表不能同时参与；避免含混地覆盖某一层。
if isfield(material,'propertyTable') && ~isempty(material.propertyTable)
    error('RocketAeroHeat:LayeredWall:Ambiguous','多层模式请把温度表放入对应层，不能同时使用顶层 propertyTable。');
end
end
