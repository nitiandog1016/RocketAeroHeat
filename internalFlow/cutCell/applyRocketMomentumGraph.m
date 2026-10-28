function image=applyRocketMomentumGraph(matrix,value)
%APPLYROCKETMOMENTUMGRAPH 完整作用原三动量A，保留跨列和法向全部原连接。
% 模块一：仅重用既有CSR，不能用法向预条件子块冒充完整A。
% 模块二：分批矩阵作用，既有代数排序不改变物理CV与原方程。
assert(strcmp(matrix.kind,'unfactored-conservative-field-block-CSR')&&matrix.dimension==3,'MomentumGraph:Kind','需要原三动量图。');
n=matrix.cellCount;assert(isequal(size(value),[n 3])&&isreal(value)&&all(isfinite(value),'all'),'MomentumGraph:Input','三动量向量无效。');
if isfield(matrix,'owners'),owners=matrix.owners;else,owners=uint32(repelem((1:n)',diff(double(matrix.rowPointer))));end
image=zeros(n,3);for first=1:4096:numel(matrix.column)
    ids=first:min(first+4095,numel(matrix.column));rows=double(owners(ids));a=rows(1);b=rows(end);
    pages=pagemtimes(matrix.values(:,:,ids),reshape(value(matrix.column(ids),:)',3,1,[]));
    for component=1:3,image(a:b,component)=image(a:b,component)+accumarray(rows-a+1,reshape(pages(component,1,:),[],1),[b-a+1 1]);end
end
assert(all(isfinite(image),'all'),'MomentumGraph:Nonfinite','原A作用非有限，拒绝使用。');
end
