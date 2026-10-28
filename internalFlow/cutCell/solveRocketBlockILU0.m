function [update,report]=solveRocketBlockILU0(factor,rhs)
%SOLVEROCKETBLOCKILU0 对固定全局五量块ILU因子作完整前后代。
% 模块一：所有控制体仅代数重排，线性右端可任意，不重筛方向或按残差早停。
% 模块二：下三角对角为I，上三角用缓存5x5逆；不声称精确求解原C或原NS。
assert(strcmp(factor.kind,'complete-five-conservation-pattern-block-ILU0-right-preconditioner'),'BlockILU:Kind','需要已验证的五量块ILU因子。');
n=factor.cellCount;assert(isequal(size(rhs),[n 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'BlockILU:RHS','五量右端非法。');
ordered=rhs(factor.order,:)';forward=zeros(5,n);solution=forward;
for row=1:n
    ids=double(factor.rowPointer(row)):double(factor.diagonalPosition(row))-1;current=ordered(:,row);
    if ~isempty(ids),product=pagemtimes(factor.values(:,:,ids),reshape(forward(:,factor.column(ids)),5,1,[]));current=current-sum(product,3);end
    forward(:,row)=current;
end
for row=n:-1:1
    ids=double(factor.diagonalPosition(row))+1:double(factor.rowPointer(row+1))-1;current=forward(:,row);
    if ~isempty(ids),product=pagemtimes(factor.values(:,:,ids),reshape(solution(:,factor.column(ids)),5,1,[]));current=current-sum(product,3);end
    solution(:,row)=factor.inverseDiagonal(:,:,row)*current;
end
update=zeros(n,5);update(factor.order,:)=solution';assert(all(isfinite(update),'all'),'BlockILU:Nonfinite','块ILU前后代非有限，不裁值。');
report=struct('status','applied-fixed-global-five-conservation-cell-block-ILU0-triangular-solves', ...
    'fixedLinearMap',true,'exactCInverseClaimed',false,'originalNSAccuracyClaimed',false, ...
    'physicalTimeAdvanced',0,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
