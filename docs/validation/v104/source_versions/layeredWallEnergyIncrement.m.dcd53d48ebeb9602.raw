function energy=layeredWallEnergyIncrement(grid,previous,current)
%LAYEREDWALLENERGYINCREMENT 逐节点材料储热变化 [J/m²]。
% 完美接触的公共节点有两个不同材料半控制体，必须分别积分 cp(T)
% 后相加；有限热阻界面本身没有热容。previous 可为统一标量初温。
if size(current,2)~=grid.nodeCount || (~isscalar(previous) && ~isequal(size(previous),size(current)))
    error('RocketAeroHeat:LayeredWall:Profile','温度剖面与多层网格不匹配。');
end
energy=zeros(size(current));
for j=1:numel(grid.layers)
    ids=grid.layerNodes{j};
    if isscalar(previous), old=previous; else, old=previous(:,ids); end
    specific=integrateWallMaterialProperty(grid.layers{j},'specificHeat',old,current(:,ids));
    energy(:,ids)=energy(:,ids)+specific.*grid.layerMass{j};
end
end
