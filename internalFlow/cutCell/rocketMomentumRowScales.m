function [scale,report]=rocketMomentumRowScales(matrix)
%ROCKETMOMENTUMROWSCALES 为完整原三动量A构造正的可逆行无穷范数尺度。
% 模块一：遍历全部原CSR系数，包含跨列与列内连接，不使用预条件近似块替代。
% 模块二：零行、尺度溢出及不可逆倒数明确拒绝，不补值或裁系数。
assert(strcmp(matrix.kind,'unfactored-conservative-field-block-CSR')&&matrix.dimension==3,'MomentumScale:Kind','需要完整原三动量A。');
n=matrix.cellCount;owners=uint32(repelem((1:n)',diff(double(matrix.rowPointer))));maximum=zeros(n,3);
for first=1:4096:numel(matrix.column)
    ids=first:min(first+4095,numel(matrix.column));rows=double(owners(ids));a=rows(1);b=rows(end);pages=max(abs(matrix.values(:,:,ids)),[],2);
    for field=1:3,maximum(a:b,field)=max(maximum(a:b,field),accumarray(rows-a+1,reshape(pages(field,1,:),[],1),[b-a+1 1],@max,0));end
end
assert(all(isfinite(maximum)&maximum>0,'all'),'MomentumScale:ZeroRow','原A有零行或非有限系数，不能构造可逆尺度。');
scale=1./maximum;assert(all(isfinite(scale)&scale>0,'all')&&all(isfinite(1./scale),'all'),'MomentumScale:Nonfinite','行尺度或逆尺度非有限。');
report=struct('method','positive-invertible-full-original-momentum-row-infinity-norm-scaling','originalRowNormRange',[min(maximum,[],'all') max(maximum,[],'all')], ...
    'rowNormSpan',max(maximum,[],'all')/min(maximum,[],'all'),'scaleRange',[min(scale,[],'all') max(scale,[],'all')], ...
    'scaledSystem','R*A*x=R*b','rightPreconditioner','P*(v./R)','originalAChanged',false,'allOriginalRowsKept',true,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
