function value=solveRocketPivotedSchur(factor,rhs)
%SOLVEROCKETPIVOTEDSCHUR 求解固定L*U≈P*S，正确应用行置换并还原全部密度/能量量。
% 模块一：接口保持两字段列序，内部可将同一CV两量相邻排列，求解后完整回映。
% 模块二：新的换行分解不放行旧2量块逆，完整C和原NS精度仍须分别验收。
assert(strcmp(factor.kind,'fixed-pivoted-scalar-ILUTP-for-density-total-energy-Schur'),'PivotedSchur:Kind','需要已固定标量主元Schur因子。');
n=factor.cellCount;assert(isequal(size(rhs),[n 2])&&isreal(rhs)&&all(isfinite(rhs),'all'),'PivotedSchur:RHS','两量右端非法。');
b=rhs(:);
hasRow=isfield(factor,'rowScale');hasColumn=isfield(factor,'columnScale');
assert(hasRow==hasColumn&&(~hasRow||isfield(factor,'scalarOrder')),'PivotedSchur:ScaleState','行列尺度及等价排列必须同时存在。');
if hasRow
    assert(isequal(size(factor.rowScale),[2*n 1])&&isequal(size(factor.columnScale),[2*n 1])&& ...
        all(isfinite(factor.rowScale)&factor.rowScale>0)&&all(isfinite(factor.columnScale)&factor.columnScale>0), ...
        'PivotedSchur:ScaleState','必须使用完整、有限、正的行列尺度。');
end
if isfield(factor,'scalarOrder')
    b=b(factor.scalarOrder);
    if isfield(factor,'rowScale'),b=factor.rowScale.*b;end
    ordered=factor.U\(factor.L\b(factor.rowPermutation));
    if isfield(factor,'columnScale'),ordered=factor.columnScale.*ordered;end
    x=zeros(2*n,1);x(factor.scalarOrder)=ordered;
else,x=factor.U\(factor.L\b(factor.rowPermutation));end
value=reshape(x,n,2);
assert(all(isfinite(value),'all'),'PivotedSchur:Nonfinite','主元Schur解非有限，不裁值。');
end
