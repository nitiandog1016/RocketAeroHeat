function policy=cutCellInteriorEulerPolicy(boundary)
%CUTCELLINTERIOREULERPOLICY 显式选择共享内面Euler通量，不在失败时偷偷切换方法。
% 模块一：缺省/显式hllc保留旧Davis-HLLC及其所有原判据。
% 模块二：exact-godunov选择精确理想气体Riemann采样，不改边界迹/黏性/墙热。
assert(isstruct(boundary)&&isscalar(boundary),'CutNS:InteriorFluxPolicy','需要标量边界结构。');
for name={'xmin','xmax','ymin','ymax','zmin','zmax','solidWall'}
    key=name{1};if isfield(boundary,key),assert(~isfield(boundary.(key),'interiorEulerFlux'), ...
        'CutNS:InteriorFluxPolicy','interiorEulerFlux须为全局选择，不能在单边或墙中被静默忽略。');end
end
policy='hllc';if isfield(boundary,'interiorEulerFlux'),policy=boundary.interiorEulerFlux;end
assert(((ischar(policy)&&isrow(policy))||(isstring(policy)&&isscalar(policy)&&~ismissing(policy)))&& ...
    ismember(char(policy),{'hllc','exact-godunov'}),'CutNS:InteriorFluxPolicy','只允许hllc或exact-godunov。');
policy=char(policy);
end
