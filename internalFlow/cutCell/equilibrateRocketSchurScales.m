function [rowScale,columnScale,report]=equilibrateRocketSchurScales(matrix,owners)
%EQUILIBRATEROCKETSCHURSCALES 构造正的可逆行列最大范数尺度，不增加或删除任何系数。
% 模块一：先R=diag(1/行∞范数)，再C=diag(1/(R*S)列∞范数)。
% 模块二：原矩阵仍为S；只对预条件分解使用R*S*C，解须回到原单位核对。
% 模块三：零行/列、溢出与非法尺度明确拒绝，不补对角或裁成有限数。
n=matrix.cellCount;nb=numel(matrix.column);assert(matrix.dimension==2&&numel(owners)==nb,'SchurScale:Input','需要完整两量图及其行编号。');
rowMax=zeros(n,2);
for first=1:4096:nb
    ids=first:min(first+4095,nb);rr=owners(ids);a=rr(1);b=rr(end);pages=abs(matrix.values(:,:,ids));maxima=max(pages,[],2);
    for field=1:2
        part=accumarray(rr-a+1,reshape(maxima(field,1,:),[],1),[b-a+1 1],@max,0);rowMax(a:b,field)=max(rowMax(a:b,field),part);
    end
end
assert(all(isfinite(rowMax),'all')&&all(rowMax>0,'all'),'SchurScale:ZeroRow','原S存在零行/非有限行，不能构造可逆尺度。');
rowScale=1./rowMax;columnMax=zeros(n,2);
assert(all(isfinite(rowScale),'all')&&all(rowScale>0,'all'),'SchurScale:Nonfinite','行尺度非正或溢出，不能继续变换。');
for first=1:4096:nb
    ids=first:min(first+4095,nb);rr=owners(ids);cc=double(matrix.column(ids));a=min(cc);b=max(cc);
    raw=abs(matrix.values(:,:,ids));pages=raw.*reshape(rowScale(rr,:)',2,1,[]);
    assert(all(isfinite(pages),'all')&&all(raw==0|pages~=0,'all'),'SchurScale:CoefficientLoss','行尺度导致非有限系数或非零系数下溢，拒绝使用。');
    maxima=max(pages,[],1);
    for field=1:2
        part=accumarray(cc-a+1,reshape(maxima(1,field,:),[],1),[b-a+1 1],@max,0);columnMax(a:b,field)=max(columnMax(a:b,field),part);
    end
end
assert(all(isfinite(columnMax),'all')&&all(columnMax>0,'all'),'SchurScale:ZeroColumn','行尺度后S存在零列/非有限列，不能补值。');
columnScale=1./columnMax;assert(all(isfinite(rowScale),'all')&&all(rowScale>0,'all')&&all(isfinite(columnScale),'all')&&all(columnScale>0,'all'), ...
    'SchurScale:Nonfinite','尺度非正或非有限，不裁值。');
report=struct('method','reversible-row-then-column-infinity-max-equilibration','rowScaleRange',[min(rowScale,[],'all') max(rowScale,[],'all')], ...
    'columnScaleRange',[min(columnScale,[],'all') max(columnScale,[],'all')],'originalRowInfinityNormRange',[min(rowMax,[],'all') max(rowMax,[],'all')], ...
    'rowScaledColumnInfinityNormRange',[min(columnMax,[],'all') max(columnMax,[],'all')], ...
    'originalSChanged',false,'coefficientsAddedOrDeleted',false,'originalResidualMustBeChecked',true,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
