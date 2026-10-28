function value=applyRocketPackedGraph(packed,update)
%APPLYROCKETPACKEDGRAPH 核对完整CSR近似矩阵C的作用，不调用原NS。
% 模块一：原CV顺序提升至代数顺序，逐行完整五量块求和，再还原物理CV。
% 模块二：该函数只适用于未因子化块图，不能把ILU因子当原矩阵。
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph'),'BlockGraph:Kind','需要未因子化的完整C块图。');
n=packed.cellCount;assert(isequal(size(update),[n 5])&&all(isfinite(update),'all'),'BlockGraph:RHS','完整块图五量输入非法。');
ordered=update(packed.order,:)';result=zeros(5,n);
for row=1:n
    ids=double(packed.rowPointer(row)):double(packed.rowPointer(row+1))-1;
    product=pagemtimes(packed.values(:,:,ids),reshape(ordered(:,packed.column(ids)),5,1,[]));result(:,row)=sum(product,3);
end
value=zeros(n,5);value(packed.order,:)=result';assert(all(isfinite(value),'all'),'BlockGraph:Nonfinite','块图作用非有限。');
end
