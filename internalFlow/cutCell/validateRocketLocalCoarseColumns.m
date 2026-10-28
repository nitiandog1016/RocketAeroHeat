function columns=validateRocketLocalCoarseColumns(columns,n,nc,np)
%VALIDATEROCKETLOCALCOARSECOLUMNS 核对冻结列身份，不按当前右端重新筛选或静默修复编号。
% 模块一：保留的父列必须递增且包括全部保守底列，附加原CV/方程对不能重复。
% 模块二：输入整数类统一为double索引；越界、缺项和重复均明确拒绝。
keys={'retainedMixedColumns';'cellIDs';'equationIDs'};
assert(isstruct(columns)&&isscalar(columns)&&isequal(sort(fieldnames(columns)),sort(keys)),'LocalCoarse:Columns','局部列配置必须明确给出父列、CV和方程身份。');
for entry=keys'
    name=entry{1};value=columns.(name);
    assert(isnumeric(value)&&isreal(value)&&isvector(value)&&~isempty(value)&&all(isfinite(value))&&all(value==floor(value)), ...
        'LocalCoarse:Columns','局部列编号必须是非空有限整数向量。');
    columns.(name)=double(value(:));
end
r=columns.retainedMixedColumns;cells=columns.cellIDs;equations=columns.equationIDs;
assert(all(r>=1&r<=nc+np)&&all(diff(r)>0)&&all(ismember((1:nc)',r))&&numel(cells)==numel(equations)&& ...
    all(cells>=1&cells<=n)&&all(equations>=1&equations<=5)&&numel(unique(cells+n*(equations-1)))==numel(cells), ...
    'LocalCoarse:Columns','冻结父列或原CV/方程索引重复、缺保守底列或越界。');
end
