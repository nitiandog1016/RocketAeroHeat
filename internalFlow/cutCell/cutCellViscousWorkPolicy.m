function policy=cutCellViscousWorkPolicy(boundary)
%CUTCELLVISCOUSWORKPOLICY 三维共享内面黏性功速度的显式全局选择。
% 模块一：缺省/legacy-wls保留全部旧数值；arithmetic只改内部面功速度。
% 模块二：单边/墙内放置该选项明确拒绝，不能静默忽略；不改牵引或Fourier。
assert(isstruct(boundary)&&isscalar(boundary),'CutNS:ViscousWorkBoundary','黏性功策略需要标量边界结构。');
for name={'xmin','xmax','ymin','ymax','zmin','zmax','solidWall'}
    key=name{1};if isfield(boundary,key),assert(~isfield(boundary.(key),'viscousWorkVelocity'), ...
        'CutNS:ViscousWorkBoundary','viscousWorkVelocity只可放在全局boundary，不能放在单边或墙内。');end
end
policy='legacy-wls';if isfield(boundary,'viscousWorkVelocity'),policy=boundary.viscousWorkVelocity;end
assert(((ischar(policy)&&isrow(policy))||(isstring(policy)&&isscalar(policy)&&~ismissing(policy)))&& ...
    ismember(char(policy),{'legacy-wls','arithmetic'}),'CutNS:ViscousWorkBoundary','黏性功策略只支持legacy-wls或arithmetic。');
policy=char(policy);
end
