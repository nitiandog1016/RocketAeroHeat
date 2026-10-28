function [image,detail]=applyRocketConsistentSchur(operator,value)
%APPLYROCKETCONSISTENTSCHUR 以同一个固定动量逆计算D-E*M_A^{-1}*B的作用。
% 模块一：与五量回代保持同一动量因子，避免将局部对角逆近似误认为该完整作用。
% 模块二：只计算已有近似C的代数作用，不改变原NS、物理网格或任何场量。
% 模块三：每次明确使用一次固定三动量前后代；不进行按右端自适应的内迭代。
assert(strcmp(operator.kind,'fixed-momentum-inverse-consistent-Schur-action'),'ConsistentSchur:Kind','需要明确绑定同一动量因子的Schur作用。');
n=operator.cellCount;assert(isequal(size(value),[n 2])&&isreal(value)&&all(isfinite(value),'all'),'ConsistentSchur:Input','需要全部CV的两量有限向量。');
bx=product(operator.B,value);momentum=solveRocketFieldILU0(operator.momentumFactor,bx);
dx=product(operator.D,value);ey=product(operator.E,momentum);image=dx-ey;
assert(all(isfinite(image),'all'),'ConsistentSchur:Nonfinite','Schur作用非有限，不裁值。');
if nargout>1,detail=struct('Bvalue',bx,'momentumValue',momentum,'Dvalue',dx,'EmomentumValue',ey, ...
    'fixedMomentumSolves',1,'blockGraphActions',3,'newActualNSEvaluations',0,'physicalEquationChanged',false);end
    function out=product(blocks,in)
        % 分批矩形块乘法：保留全部原连接，仅在已有代数排序内累加。
        dimension=size(blocks,1);d=size(in,2);out=zeros(n,dimension);nb=numel(operator.column);
        for first=1:4096:nb
            ids=first:min(first+4095,nb);rows=double(operator.owners(ids));a=rows(1);b=rows(end);
            pages=pagemtimes(blocks(:,:,ids),reshape(in(operator.column(ids),:)',d,1,[]));
            for component=1:dimension,out(a:b,component)=out(a:b,component)+accumarray(rows-a+1,reshape(pages(component,1,:),[],1),[b-a+1 1]);end
        end
    end
end
