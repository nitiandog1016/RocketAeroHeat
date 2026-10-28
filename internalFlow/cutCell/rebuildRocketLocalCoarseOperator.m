function coarse=rebuildRocketLocalCoarseOperator(cache,space,primitiveSpace,selected,columns,maximumImageBytes)
%REBUILDROCKETLOCALCOARSEOPERATOR 新Q只按冻结列身份重算当前像，原保守底空间完整保留。
% 模块一：仅构造保留父列和指定原CV/方程，不重筛右端或构造已释放父列。
% 模块二：每列重新归一化，完整Gram正定/条件与原预算失败时明确拒绝。
n=size(cache.Q,1);nc=5*space.groupCount;np=numel(selected);selected=selected(:);
columns=validateRocketLocalCoarseColumns(columns,n,nc,np);nr=numel(columns.retainedMixedColumns);nl=numel(columns.cellIDs);d=nr+nl;
assert(numel(unique(selected))==np&&all(selected>=1&selected<=5*primitiveSpace.groupCount&selected==floor(selected)), ...
    'LocalCoarse:Columns','父混合原始变量列非法。');
mixed=struct('space',struct('groupCount',space.groupCount,'prolongation',space.prolongation), ...
    'primitiveSpace',rmfield(primitiveSpace,'geometry'),'baseQ',cache.Q,'gas',cache.gas,'scale',cache.scale, ...
    'selectedPrimitiveColumns',selected,'columnLiftScale',ones(nc+np,1));
nf=size(cache.geometry.inputGeometry.F,1);
coarse=struct('kind','local-conservative-enriched-mixed-fine-L2-minimum-residual-coarse-correction', ...
    'space',space,'baseQ',cache.Q,'gas',cache.gas,'boundary',cache.boundary,'scale',cache.scale, ...
    'sourceIdentity',rocketSelectedDirectionSourceIdentity,'liftSourceIdentity',rocketLocalCoarseLiftIdentity,'mixedLift',mixed, ...
    'retainedMixedColumns',columns.retainedMixedColumns,'baseConservativeColumns',nc, ...
    'localCellIDs',columns.cellIDs,'localEquationIDs',columns.equationIDs,'localColumnLiftScale',ones(nl,1), ...
    'localOriginalFaceIDs',mod(columns.cellIDs-1,nf)+1,'localLayerIDs',floor((columns.cellIDs-1)/nf)+1);
imageColumns=cell(1,d);liftScale=zeros(d,1);count=0;timer=tic;
for col=1:d
    volumeScale=1;
    if col<=nr
        identity=columns.retainedMixedColumns(col);
        if identity<=nc,g=mod(identity-1,space.groupCount)+1;volumeScale=1/sqrt(space.groupVolume(g));
        else,g=mod(selected(identity-nc)-1,primitiveSpace.groupCount)+1;volumeScale=1/sqrt(primitiveSpace.groupVolume(g));end
    end
    weights=zeros(d,1);weights(col)=volumeScale;direction=liftRocketLocalCoarse(coarse,weights);
    value=-applyCutCellSelectedDirection(cache,direction);magnitude=norm(value(:));assert(isfinite(magnitude)&&magnitude>0,'LocalCoarse:NullImage','新Q的局部列像无效。');
    imageColumns{col}=sparse(value(:)/magnitude);liftScale(col)=volumeScale/magnitude;count=count+nnz(imageColumns{col});
    assert(16*count+8*(d+1)<=maximumImageBytes,'LocalCoarse:Budget','新Q局部像超过原预算，不裁细像或删物理CV。');
end
image=horzcat(imageColumns{:});clear imageColumns;buffers=whos('image');assert(buffers.bytes<=maximumImageBytes,'LocalCoarse:Budget','新Q实际像超界。');
G=full(image'*image);G=.5*(G+G');condition=rcond(G);[upper,flag]=chol(G);
assert(isfinite(condition)&&condition>1e-12&&flag==0,'LocalCoarse:Factor','新Q局部Gram非正定或病态，不补对角。');
coarse.mixedLift.columnLiftScale(columns.retainedMixedColumns)=liftScale(1:nr);coarse.localColumnLiftScale=liftScale(nr+1:end);coarse.columnLiftScale=liftScale;
coarse.image=image;coarse.normalMatrix=G;coarse.upper=upper;coarse.imageBytes=buffers.bytes;coarse.imageNonzeros=nnz(image);
coarse.maximumImageBytes=maximumImageBytes;coarse.reciprocalConditionEstimate=condition;coarse.constructionModelApplications=d;coarse.copiedSourceImageColumns=0;
coarse.buildSeconds=toc(timer);coarse.actualOriginalNSEvaluations=0;coarse.fineJacobianMatrixAllocated=false;coarse.newProductionSolvedAreaM2=0;
coarse.selectionStopReason='fixed-local-and-retained-mixed-columns-rebuilt-at-current-Q';
end
