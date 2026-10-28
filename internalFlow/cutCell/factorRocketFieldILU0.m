function factor=factorRocketFieldILU0(factor,observation)
%FACTORROCKETFIELDILU0 对三动量或密度能量子块作固定模式ILU，保留全部原CV。
% 模块一：仅处理2×2/3×3完整对角块，块内不假定标量三角性。
% 模块二：不补对角；奇异或原1e-8局部逆检查失败明确拒绝。
% 模块三：显式观察器只在原拒绝前保存主元/消元账及不可用于求解的部分因子。
if nargin<2,observation=struct([]);end
observe=~isempty(observation);
if observe
    assert(isstruct(observation)&&isscalar(observation)&&all(isfield(observation, ...
        {'failureFile','sourceBinding','assemblyFile','cellOrder','sourceFaces','fieldGroup'})), ...
        'FieldILU:Observation','失败观察器缺少明确路径/来源/CV映射。');
    assert(ischar(observation.failureFile)&&exist(observation.failureFile,'file')~=2&& ...
        exist(observation.assemblyFile,'file')==2,'FieldILU:Observation','必须使用新失败文件和已保存组装检查点。');
end
assert(strcmp(factor.kind,'unfactored-conservative-field-block-CSR')&&ismember(factor.dimension,[2 3]), ...
    'FieldILU:Kind','需要未因子化的2或3量子块图。');
n=factor.cellCount;d=factor.dimension;inverse=zeros(d,d,n);position=zeros(n,1,'uint32');minimumRcond=Inf;maximumDefect=0;retained=0;discarded=0;timer=tic;
for row=1:n
    ids=double(factor.rowPointer(row)):double(factor.rowPointer(row+1))-1;columns=double(factor.column(ids));position(columns)=uint32(ids);diagonal=double(factor.diagonalPosition(row));
    if observe,rowInput=factor.values(:,:,ids);end
    for entry=double(factor.rowPointer(row)):diagonal-1
        pivot=double(factor.column(entry));mult=factor.values(:,:,entry)*inverse(:,:,pivot);factor.values(:,:,entry)=mult;
        up=double(factor.diagonalPosition(pivot))+1:double(factor.rowPointer(pivot+1))-1;where=position(factor.column(up));keep=where>0;
        retained=retained+nnz(keep);discarded=discarded+nnz(~keep);
        if any(keep),factor.values(:,:,double(where(keep)))=factor.values(:,:,double(where(keep)))-pagemtimes(mult,factor.values(:,:,up(keep)));end
    end
    D=factor.values(:,:,diagonal);condition=rcond(D);
    if observe&&~(isfinite(condition)&&condition>0)
        saveFailure(factor,inverse,rowInput,row,D,nan(d),condition,NaN,observation,'FieldILU:Singular');
    end
    assert(isfinite(condition)&&condition>0,'FieldILU:Singular','%d量子块第%d主元奇异，不补对角。',d,row);
    inv=D\eye(d);defect=max(abs(D*inv-eye(d)),[],'all');
    if observe&&~(all(isfinite(inv),'all')&&defect<1e-8)
        saveFailure(factor,inverse,rowInput,row,D,inv,condition,defect,observation,'FieldILU:FactorDefect');
    end
    assert(all(isfinite(inv),'all')&&defect<1e-8,'FieldILU:FactorDefect','%d量子块第%d局部逆未过原1e-8。',d,row);
    inverse(:,:,row)=inv;minimumRcond=min(minimumRcond,condition);maximumDefect=max(maximumDefect,defect);position(columns)=0;
    if mod(row,50000)==0,fprintf('FIELD_ILU d=%d row=%d/%d\n',d,row,n);end
end
factor.inverseDiagonal=inverse;factor.kind='factored-conservative-field-block-ILU0';
factor.factorizationReport=struct('dimension',d,'minimumLocalPivotRcond',minimumRcond,'maximumLocalInverseDefect',maximumDefect, ...
    'retainedProducts',retained,'discardedFillProducts',discarded,'factorSeconds',toc(timer),'diagonalCompensationApplied',false,'allOriginalCellsKept',true);
end
function saveFailure(partial,cachedInverse,rowInput,row,D,inverseCandidate,condition,defect,observation,errorIdentifier)
% 仅重放被拒绝行的已有消元，保存原始行、每个保留乘积与真实被拒主元。
n=partial.cellCount;d=partial.dimension;first=double(partial.rowPointer(row));last=double(partial.rowPointer(row+1))-1;
diagonal=double(partial.diagonalPosition(row));columns=double(partial.column(first:last));position=zeros(n,1,'uint32');position(columns)=uint32(first:last);
work=rowInput;ledger={};formation=norm(rowInput(:,:,diagonal-first+1),'fro');
for entry=first:diagonal-1
    pivot=double(partial.column(entry));before=work(:,:,entry-first+1);mult=before*cachedInverse(:,:,pivot);work(:,:,entry-first+1)=mult;
    upper=double(partial.diagonalPosition(pivot))+1:double(partial.rowPointer(pivot+1))-1;where=position(partial.column(upper));keep=where>0;
    kept=upper(keep);targets=double(where(keep));upperBlocks=partial.values(:,:,kept);products=pagemtimes(mult,upperBlocks);beforeD=work(:,:,diagonal-first+1);
    if ~isempty(targets),work(:,:,targets-first+1)=work(:,:,targets-first+1)-products;end
    selected=find(targets==diagonal);if ~isempty(selected),formation=formation+norm(products(:,:,selected),'fro');end
    ledger{end+1}=struct('pivotAlgebraicRow',pivot,'lowerEntry',entry,'lowerBlockBefore',before, ...
        'cachedPivotInverse',cachedInverse(:,:,pivot),'multiplier',mult,'upperEntries',kept,'targetEntries',targets, ...
        'upperBlocks',upperBlocks,'retainedProducts',products,'discardedColumns',partial.column(upper(~keep)), ...
        'diagonalBefore',beforeD,'diagonalAfter',work(:,:,diagonal-first+1)); %#ok<AGROW>
end
actual=partial.values(:,:,first:last);originalCV=observation.cellOrder(row);sourceFace=mod(originalCV-1,observation.sourceFaces)+1;layer=floor((originalCV-1)/observation.sourceFaces)+1;
partial.kind='partial-field-ILU-before-rejected-inverse-commit';
failure=struct('status','rejected-by-unchanged-original-local-inverse-gate','errorIdentifier',errorIdentifier, ...
    'dimension',d,'fieldGroup',observation.fieldGroup,'algebraicRow',row,'lastCommittedInverseRow',row-1, ...
    'originalControlVolume',originalCV,'originalFace',sourceFace,'layer',layer,'pivotD',D,'inverseCandidate',inverseCandidate, ...
    'inverseProduct',D*inverseCandidate,'maximumInverseDefect',defect,'pivotRcond',condition,'originalInverseDefectThreshold',1e-8, ...
    'rowColumns',columns,'rowInput',rowInput,'retainedEliminationLedger',{ledger},'replayedRow',work,'actualPartialRow',actual, ...
    'ledgerExactlyReproducesPartialRow',isequaln(work,actual),'ledgerRelativeDifference',norm(work(:)-actual(:))/max(1,norm(actual(:))), ...
    'diagonalFormationNormBound',formation,'diagonalCancellationRatio',formation/max(realmin,norm(D,'fro')), ...
    'partialFactor',partial,'cachedInverseDiagonal',cachedInverse,'sourceBinding',observation.sourceBinding, ...
    'assemblyFile',observation.assemblyFile,'assemblyFileSHA256',rocketStreamingFileSHA256(observation.assemblyFile), ...
    'preconditionerUsable',false,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
save(observation.failureFile,'failure','-v7');
summary=rmfield(failure,{'partialFactor','cachedInverseDiagonal','retainedEliminationLedger','rowInput','replayedRow','actualPartialRow'});
summary.failureFile=observation.failureFile;summary.failureFileSHA256=rocketStreamingFileSHA256(observation.failureFile);
summary.retainedEliminationSteps=numel(ledger);file=[observation.failureFile '.json'];assert(exist(file,'file')~=2,'FieldILU:Observation','不得覆盖失败摘要。');
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(summary,'ConvertInfAndNaN',true));
end
