function coarse=rebuildRocketMixedCoarseOperator(cache,space,primitiveSpace,selected,maximumImageBytes)
%REBUILDROCKETMIXEDCOARSEOPERATOR 接受新Q后按固定列身份重算混合像，不复用旧数值或重筛选右端。
% 模块一：保守底空间所有列与所选原始变量列各执行一次当前完整模型。
% 模块二：新像重新归一化/分解，条件或预算失败明确拒绝，不裁值/补对角。
selected=selected(:);nc=5*space.groupCount;d=nc+numel(selected);n=size(cache.Q,1);
assert(numel(unique(selected))==numel(selected)&&all(selected>=1&selected<=5*primitiveSpace.groupCount&selected==floor(selected)),'MixedCoarse:Columns','混合列身份非法。');
coarse=struct('kind','mixed-conservative-primitive-fine-L2-minimum-residual-coarse-correction','space',space, ...
    'primitiveSpace',rmfield(primitiveSpace,'geometry'),'baseQ',cache.Q,'gas',cache.gas,'boundary',cache.boundary,'scale',cache.scale, ...
    'sourceIdentity',rocketSelectedDirectionSourceIdentity,'liftSourceIdentity',rocketMixedCoarseLiftIdentity,'selectedPrimitiveColumns',selected, ...
    'baseConservativeColumns',nc,'columnLiftScale',ones(d,1));
columns=cell(1,d);liftScale=zeros(d,1);count=0;timer=tic;
for col=1:d
    if col<=nc,g=mod(col-1,space.groupCount)+1;vscale=1/sqrt(space.groupVolume(g));
    else,g=mod(selected(col-nc)-1,primitiveSpace.groupCount)+1;vscale=1/sqrt(primitiveSpace.groupVolume(g));end
    coefficients=zeros(d,1);coefficients(col)=vscale;direction=liftRocketMixedCoarse(coarse,coefficients);
    value=-applyCutCellSelectedDirection(cache,direction);magnitude=norm(value(:));assert(isfinite(magnitude)&&magnitude>0,'MixedCoarse:NullImage','当前混合方向像无效。');
    columns{col}=sparse(value(:)/magnitude);liftScale(col)=vscale/magnitude;count=count+nnz(columns{col});
    assert(16*count+8*(d+1)<=maximumImageBytes,'MixedCoarse:Budget','新Q混合像超过原预算，不截非零。');
end
image=horzcat(columns{:});clear columns;buffers=whos('image');assert(buffers.bytes<=maximumImageBytes,'MixedCoarse:Budget','新Q实际混合像超过预算。');
G=full(image'*image);G=.5*(G+G');condition=rcond(G);assert(isfinite(condition)&&condition>1e-12,'MixedCoarse:Condition','新Q混合Gram病态，不补对角。');
[upper,flag]=chol(G);assert(flag==0,'MixedCoarse:Factor','新Q混合Gram非正定。');
coarse.columnLiftScale=liftScale;coarse.image=image;coarse.normalMatrix=G;coarse.upper=upper;coarse.imageBytes=buffers.bytes;coarse.imageNonzeros=nnz(image);
coarse.maximumImageBytes=maximumImageBytes;coarse.reciprocalConditionEstimate=condition;coarse.constructionModelApplications=d;coarse.copiedSourceImageColumns=0;
coarse.buildSeconds=toc(timer);coarse.actualOriginalNSEvaluations=0;coarse.fineJacobianMatrixAllocated=false;coarse.newProductionSolvedAreaM2=0;
coarse.selectionStopReason='fixed-selected-basis-rebuilt-at-current-Q';coarse.selectionHistory=[];
end
