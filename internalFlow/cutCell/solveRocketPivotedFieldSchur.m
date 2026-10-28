function [value,report]=solveRocketPivotedFieldSchur(pre,rhs)
%SOLVEROCKETPIVOTEDFIELDSCHUR 用三动量块解与行主元密度/能量Schur组成固定五量映射。
% 模块一：完整B/E与两次动量回代保持，全部原CV/方程只作可逆代数排列。
% 模块二：稀疏丢弃仅影响预条件；不修改原C/NS，也不将此结果作为气动热解。
assert(strcmp(pre.kind,'fixed-five-conservation-with-pivoted-scalar-Schur'),'PivotedFieldSchur:Kind','需要固定五量主元Schur预条件器。');
n=pre.cellCount;assert(isequal(size(rhs),[n 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'PivotedFieldSchur:RHS','五量右端非法。');
r=rhs(pre.order,:);ym=solveRocketFieldILU0(pre.momentumFactor,r(:,2:4));rs=r(:,[1 5])-coupling(pre.E,ym,2);
ys=solveRocketPivotedSchur(pre.schurFactor,rs);ym=ym-solveRocketFieldILU0(pre.momentumFactor,coupling(pre.B,ys,3));
sorted=zeros(n,5);sorted(:,2:4)=ym;sorted(:,[1 5])=ys;value=zeros(n,5);value(pre.order,:)=sorted;
assert(all(isfinite(value),'all'),'PivotedFieldSchur:Nonfinite','五量主元Schur解非有限。');
report=struct('method','fixed-three-momentum-ILU-and-row-pivoted-density-energy-Schur','fixedLinearMap',true,'oldBlockInverseClaimedValid',false, ...
    'newActualNSEvaluations',0,'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    function out=coupling(blocks,in,dimension)
        out=zeros(n,dimension);d=size(in,2);nb=numel(pre.column);
        for first=1:4096:nb
            ids=first:min(first+4095,nb);rows=double(pre.owners(ids));a=rows(1);b=rows(end);products=pagemtimes(blocks(:,:,ids),reshape(in(pre.column(ids),:)',d,1,[]));
            for component=1:dimension,out(a:b,component)=out(a:b,component)+accumarray(rows-a+1,reshape(products(component,1,:),[],1),[b-a+1 1]);end
        end
    end
end
