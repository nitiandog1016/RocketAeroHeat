function value=solveRocketFieldILU0(factor,rhs)
%SOLVEROCKETFIELDILU0 在已固定代数顺序中前后代，完整对角块由缓存逆求解。
% 模块一：子块保持全部CV，不作残差早停或按右端裁剪；因此是固定线性映射。
% 模块二：上三角是块上三角，块内所有分量必须参与求解。
assert(strcmp(factor.kind,'factored-conservative-field-block-ILU0'),'FieldILU:Kind','需要已因子化子块。');
n=factor.cellCount;d=factor.dimension;assert(isequal(size(rhs),[n d])&&isreal(rhs)&&all(isfinite(rhs),'all'),'FieldILU:RHS','子块右端非法。');
b=rhs';forward=zeros(d,n);solution=forward;
for row=1:n
    ids=double(factor.rowPointer(row)):double(factor.diagonalPosition(row))-1;current=b(:,row);
    if ~isempty(ids),current=current-sum(pagemtimes(factor.values(:,:,ids),reshape(forward(:,factor.column(ids)),d,1,[])),3);end
    forward(:,row)=current;
end
for row=n:-1:1
    ids=double(factor.diagonalPosition(row))+1:double(factor.rowPointer(row+1))-1;current=forward(:,row);
    if ~isempty(ids),current=current-sum(pagemtimes(factor.values(:,:,ids),reshape(solution(:,factor.column(ids)),d,1,[])),3);end
    solution(:,row)=factor.inverseDiagonal(:,:,row)*current;
end
value=solution';assert(all(isfinite(value),'all'),'FieldILU:Nonfinite','子块前后代非有限，不裁值。');
end
