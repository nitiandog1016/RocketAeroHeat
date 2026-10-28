function value=applyRocketPackedGraphTranspose(packed,update)
%APPLYROCKETPACKEDGRAPHTRANSPOSE 对未因子化完整五量CSR作精确转置作用。
% 模块一：按原CV/代数排列变换，保留每个原块，不对称化或删去跨组连接。
% 模块二：分批计算每个原块的转置贡献，再按原列累加，避免形成全局稠密矩阵。
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph'),'BlockGraph:Kind','需要未因子化的完整五量块图。');
n=packed.cellCount;assert(isequal(size(update),[n 5])&&isreal(update)&&all(isfinite(update),'all'),'BlockGraph:RHS','完整图转置五量输入非法。');
ordered=update(packed.order,:);result=zeros(n,5);owners=repelem((1:n)',diff(double(packed.rowPointer)));
for first=1:32768:numel(packed.column)
    ids=first:min(numel(packed.column),first+32767);b=reshape(ordered(owners(ids),:)',5,1,[]);
    product=pagemtimes(permute(packed.values(:,:,ids),[2 1 3]),b);product=reshape(product,5,[])';
    column=double(packed.column(ids));
    for field=1:5,result(:,field)=result(:,field)+accumarray(column(:),product(:,field),[n 1]);end
end
value=zeros(n,5);value(packed.order,:)=result;
assert(all(isfinite(value),'all'),'BlockGraph:Nonfinite','完整块图转置作用非有限。');
end
