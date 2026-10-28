function [value,report]=solveRocketFieldSchurPreconditioner(pre,rhs)
%SOLVEROCKETFIELDSCHURPRECONDITIONER 固定五量Schur分块前后代，不求新物理状态。
% 模块一：先三动量、再密度/能量Schur、最后回代三动量，保留全部交叉耦合B/E。
% 模块二：S采用局部动量逆构造而非精确全局Schur，完整C/原NS精度仍需单独验证。
assert(strcmp(pre.kind,'fixed-five-conservation-momentum-density-energy-approximate-Schur'),'FieldSchur:Kind','需要已固定五量Schur预条件器。');
n=pre.cellCount;assert(isequal(size(rhs),[n 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'FieldSchur:RHS','五量右端非法。');
r=rhs(pre.order,:);ym=solveRocketFieldILU0(pre.momentumFactor,r(:,2:4));rs=r(:,[1 5])-coupling(pre.E,ym,2);
ys=solveRocketFieldILU0(pre.schurFactor,rs);ym=ym-solveRocketFieldILU0(pre.momentumFactor,coupling(pre.B,ys,3));
sorted=zeros(n,5);sorted(:,2:4)=ym;sorted(:,[1 5])=ys;value=zeros(n,5);value(pre.order,:)=sorted;
assert(all(isfinite(value),'all'),'FieldSchur:Nonfinite','Schur分块解非有限，不裁值。');
report=struct('status','applied-fixed-five-conservation-approximate-Schur-block-factorization', ...
    'momentumTriangularSolves',2,'densityEnergyTriangularSolves',1,'fixedLinearMap',true,'exactCInverseClaimed',false, ...
    'newActualNSEvaluations',0,'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    function out=coupling(blocks,in,dimension)
        out=zeros(n,dimension);d=size(in,2);nb=numel(pre.column);
        for first=1:4096:nb
            ids=first:min(first+4095,nb);rows=double(pre.owners(ids));a=rows(1);b=rows(end);
            products=pagemtimes(blocks(:,:,ids),reshape(in(pre.column(ids),:)',d,1,[]));
            for component=1:dimension
                out(a:b,component)=out(a:b,component)+accumarray(rows-a+1,reshape(products(component,1,:),[],1),[b-a+1 1]);
            end
        end
    end
end
