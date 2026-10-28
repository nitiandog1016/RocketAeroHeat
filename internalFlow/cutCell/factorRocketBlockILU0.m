function factor=factorRocketBlockILU0(factor)
%FACTORROCKETBLOCKILU0 在完整五量控制体连接上做块ILU(0)，保留全部原矩阵块。
% 模块一：仅舍弃消元产生的模式外填充，不删除原CV/方程或原矩阵连接。
% 模块二：5x5主元内部使用数值求逆，奇异或原1e-8逆缺陷门槛失败直接拒绝。
% 模块三：不补对角、不改变C或原NS；因子只用于固定线性右预条件。
assert(strcmp(factor.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph'),'BlockILU:Kind','需要未因子化的完整块图。');
n=factor.cellCount;inverseDiagonal=zeros(5,5,n);position=zeros(n,1,'uint32');identity=eye(5);
minimumRcond=Inf;maximumDefect=0;retainedProducts=0;discardedFillProducts=0;timer=tic;
for row=1:n
    ids=double(factor.rowPointer(row)):double(factor.rowPointer(row+1))-1;
    column=double(factor.column(ids));position(column)=uint32(ids);diagonal=double(factor.diagonalPosition(row));
    for entry=double(factor.rowPointer(row)):diagonal-1
        pivot=double(factor.column(entry));assert(pivot<row,'BlockILU:Pattern','下带排序无效。');
        multiplier=factor.values(:,:,entry)*inverseDiagonal(:,:,pivot);factor.values(:,:,entry)=multiplier;
        upper=double(factor.diagonalPosition(pivot))+1:double(factor.rowPointer(pivot+1))-1;
        where=position(factor.column(upper));keep=where>0;retainedProducts=retainedProducts+nnz(keep);discardedFillProducts=discardedFillProducts+nnz(~keep);
        if any(keep)
            target=double(where(keep));change=pagemtimes(multiplier,factor.values(:,:,upper(keep)));
            factor.values(:,:,target)=factor.values(:,:,target)-change;
        end
    end
    D=factor.values(:,:,diagonal);condition=rcond(D);
    assert(isfinite(condition)&&condition>0,'BlockILU:Singular','第%d代数主元（原CV%d）奇异，不补对角。',row,factor.order(row));
    inverse=D\identity;defect=max(abs(D*inverse-identity),[],'all');
    assert(all(isfinite(inverse),'all')&&defect<1e-8,'BlockILU:FactorDefect','第%d代数主元逆缺陷未通过原门槛。',row);
    inverseDiagonal(:,:,row)=inverse;minimumRcond=min(minimumRcond,condition);maximumDefect=max(maximumDefect,defect);position(column)=0;
    if mod(row,10000)==0,fprintf('BLOCK_ILU row=%d/%d retained=%d minBlockRcond=%.6g\n',row,n,retainedProducts,minimumRcond);end
end
factor.inverseDiagonal=inverseDiagonal;factor.kind='complete-five-conservation-pattern-block-ILU0-right-preconditioner';
factor.factorizationReport=struct('method','five-by-five-cell-block-ILU-zero-fill','retainedBlockProducts',retainedProducts, ...
    'discardedOutOfPatternFillProducts',discardedFillProducts,'minimumLocalPivotRcond',minimumRcond, ...
    'maximumLocalInverseDefect',maximumDefect,'factorSeconds',toc(timer),'diagonalCompensationApplied',false, ...
    'allOriginalCellsAndFiveEquationsKept',true,'originalGraphConnectionsDropped',false,'fixedLinearMap',true, ...
    'physicalEquationChanged',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
