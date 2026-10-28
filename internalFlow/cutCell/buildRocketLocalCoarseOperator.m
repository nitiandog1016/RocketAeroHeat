function coarse=buildRocketLocalCoarseOperator(parent,cache,remaining,cellIDs,options)
%BUILDROCKETLOCALCOARSEOPERATOR 用原CV局部保守方向补充混合空间，最终像遵守原内存预算。
% 模块一：每个候选是一个原CV/保守量的单位方向，完整当前空间作用不裁非零。
% 模块二：必要时仅释放低投影收益的整列原始变量方向，保留全部保守底空间列。
% 模块三：Schur收益筛选和条件检查只在构造时执行；后续右端使用冻结列。
if nargin<5,options=struct;end
defaults=struct('maximumImageBytes',256*1024^2,'maximumAdditionalColumns',320, ...
    'minimumOrthogonalImageSquaredNorm',1e-8,'minimumRelativeGain',1e-10);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'LocalCoarse:Options','未知局部粗层选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumAdditionalColumns,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
for key={'maximumImageBytes','minimumOrthogonalImageSquaredNorm','minimumRelativeGain'},validateattributes(options.(key{1}),{'numeric'},{'real','finite','scalar','positive'});end
assert(strcmp(parent.kind,'mixed-conservative-primitive-fine-L2-minimum-residual-coarse-correction')&& ...
    isequaln(parent.baseQ,cache.Q)&&isequaln(parent.space.geometry,cache.geometry)&&isequaln(parent.gas,cache.gas)&& ...
    isequaln(parent.boundary,cache.boundary)&&isequaln(parent.scale,cache.scale)&& ...
    isequaln(parent.sourceIdentity,rocketSelectedDirectionSourceIdentity)&&isequaln(parent.liftSourceIdentity,rocketMixedCoarseLiftIdentity), ...
    'LocalCoarse:Parent','父像基态、空间或来源与当前模型不符。');
assert(isequal(size(remaining),size(cache.Q))&&all(isfinite(remaining),'all'),'LocalCoarse:RHS','筛选剩余非法。');
n=size(cache.Q,1);validateattributes(cellIDs,{'numeric'},{'real','finite','vector','integer','positive','<=',n});cellIDs=double(cellIDs(:));
assert(numel(unique(cellIDs))==numel(cellIDs),'LocalCoarse:Cells','候选原CV不能重复。');
timer=tic;candidateCells=repelem(cellIDs,5);candidateEquations=repmat((1:5)',numel(cellIDs),1);np=numel(candidateCells);
columns=cell(1,np);liftScale=zeros(np,1);
for col=1:np
    direction=zeros(size(cache.Q));direction(candidateCells(col),candidateEquations(col))=1;
    value=-applyCutCellSelectedDirection(cache,direction);magnitude=norm(value(:));
    assert(isfinite(magnitude)&&magnitude>0,'LocalCoarse:NullImage','局部单位方向没有有效的完整像。');
    columns{col}=sparse(value(:)/magnitude);liftScale(col)=1/magnitude;
    if mod(col,25)==0||col==np,fprintf('LOCAL_COARSE_IMAGE columns=%d/%d\n',col,np);end
end
L=horzcat(columns{:});clear columns;candidateNNZ=full(sum(spones(L),1))';
S=parent.image;Gparent=full(S'*S);Gparent=.5*(Gparent+Gparent');
assert(norm(Gparent-parent.normalMatrix,'fro')<1e-10*max(1,norm(Gparent,'fro'))&&rcond(Gparent)>1e-12,'LocalCoarse:Parent','父像小Gram校核失败。');
retained=(1:size(S,2))';removed=zeros(0,1);pruneHistory=zeros(0,4);vector=remaining(:);sourceNNZ=full(sum(spones(S),1))';
% 以全部候选可容纳为上界释放整列；不截细像、不把256MiB称进程峰值内存。
while 16*(sum(sourceNNZ(retained))+sum(candidateNNZ))+8*(numel(retained)+np+1)>options.maximumImageBytes
    G=Gparent(retained,retained);inverse=G\eye(numel(retained));weights=G\(S(:,retained)'*vector);
    loss=weights.^2./diag(inverse);available=find(retained>parent.baseConservativeColumns);
    assert(~isempty(available),'LocalCoarse:Budget','全部局部候选与保守底空间超过预算，不删除物理单元或裁细像。');
    [~,which]=min(loss(available)./(16*sourceNNZ(retained(available))+8));position=available(which);old=retained(position);
    pruneHistory(end+1,:)=[old loss(position) 16*sourceNNZ(old)+8 numel(retained)-1]; %#ok<AGROW>
    removed(end+1,1)=old;retained(position)=[]; %#ok<AGROW>
end
B=S(:,retained);G0=Gparent(retained,retained);[upper,flag]=chol(G0);assert(flag==0&&rcond(G0)>1e-12,'LocalCoarse:Base','保留父空间条件失败。');
cross=full(B'*L);GP=full(L'*L);Y=upper\(upper'\cross);H=GP-cross'*Y;H=.5*(H+H');
baseWeights=upper\(upper'\(B'*vector));baseRemaining=vector-B*baseWeights;correlation0=L'*baseRemaining;
selected=zeros(0,1);excluded=false(np,1);reasons=repmat({''},np,1);history=zeros(0,7);nnzCount=nnz(B);G=G0;stop='maximum-additional-columns';
for attempt=1:np
    if numel(selected)>=options.maximumAdditionalColumns,break;end
    if isempty(selected),diagonal=diag(H);correlation=correlation0;
    else
        small=H(selected,selected);diagonal=diag(H)-sum(H(:,selected).*(small\H(selected,:))',2);
        correlation=correlation0-H(:,selected)*(small\correlation0(selected));
    end
    valid=~excluded;valid(selected)=false;independent=diagonal>options.minimumOrthogonalImageSquaredNorm;
    newly=valid&~independent;excluded(newly)=true;if any(newly),reasons(newly)={'dependent-at-selection-tolerance'};end
    score=-inf(np,1);valid=valid&independent;score(valid)=correlation(valid).^2./diagonal(valid);
    [gain,j]=max(score);
    if ~isfinite(gain),stop='no-independent-candidate';break;end
    if gain<=options.minimumRelativeGain*max(1,sum(vector.^2)),stop='additional-projection-gain-below-construction-threshold';break;end
    trial=[selected;j];trialG=[G0 cross(:,trial);cross(:,trial)' GP(trial,trial)];trialG=.5*(trialG+trialG');
    condition=rcond(trialG);[factor,flag]=chol(trialG);
    if flag~=0||~isfinite(condition)||condition<=1e-12,excluded(j)=true;reasons{j}='full-Gram-conditioning-or-positive-definiteness';continue;end
    selected=trial;G=trialG;upper=factor;nnzCount=nnzCount+candidateNNZ(j);
    weights=upper\(upper'\[B'*vector;L(:,selected)'*vector]);predicted=vector-B*weights(1:numel(retained))-L(:,selected)*weights(numel(retained)+1:end);
    history(end+1,:)=[numel(selected) j diagonal(j) gain condition 16*nnzCount+8*(numel(retained)+numel(selected)+1) norm(predicted)/max(1,norm(vector))]; %#ok<AGROW>
end
image=[B L(:,selected)];buffers=whos('image');assert(buffers.bytes<=options.maximumImageBytes,'LocalCoarse:Budget','最终局部像超过原预算。');
assert(norm(full(image'*image)-G,'fro')<1e-10*max(1,norm(G,'fro')),'LocalCoarse:Gram','最终像与Gram不符。');
mixedLift=struct('space',struct('groupCount',parent.space.groupCount,'prolongation',parent.space.prolongation), ...
    'primitiveSpace',parent.primitiveSpace,'baseQ',parent.baseQ,'gas',parent.gas,'scale',parent.scale, ...
    'selectedPrimitiveColumns',parent.selectedPrimitiveColumns,'columnLiftScale',parent.columnLiftScale);
nf=size(cache.geometry.inputGeometry.F,1);
coarse=struct('kind','local-conservative-enriched-mixed-fine-L2-minimum-residual-coarse-correction', ...
    'space',parent.space,'baseQ',parent.baseQ,'gas',parent.gas,'boundary',parent.boundary,'scale',parent.scale, ...
    'sourceIdentity',parent.sourceIdentity,'liftSourceIdentity',rocketLocalCoarseLiftIdentity,'mixedLift',mixedLift, ...
    'retainedMixedColumns',retained,'removedMixedColumns',removed,'baseConservativeColumns',parent.baseConservativeColumns, ...
    'localCellIDs',candidateCells(selected),'localEquationIDs',candidateEquations(selected),'localColumnLiftScale',liftScale(selected), ...
    'localOriginalFaceIDs',mod(candidateCells(selected)-1,nf)+1,'localLayerIDs',floor((candidateCells(selected)-1)/nf)+1, ...
    'columnLiftScale',[parent.columnLiftScale(retained);liftScale(selected)],'image',image,'normalMatrix',G,'upper',upper, ...
    'imageBytes',buffers.bytes,'imageNonzeros',nnz(image),'maximumImageBytes',options.maximumImageBytes,'reciprocalConditionEstimate',rcond(G), ...
    'selectionOptions',options,'selectionRemaining',remaining,'pruneHistory',pruneHistory,'selectionHistory',history, ...
    'selectionStopReason',stop,'excludedCandidateReasons',{reasons},'candidateCellIDs',candidateCells,'candidateEquationIDs',candidateEquations, ...
    'candidateImage',L,'candidateColumnLiftScale',liftScale,'selectedCandidateColumns',selected, ...
    'constructionModelApplications',np,'copiedSourceImageColumns',numel(retained),'buildSeconds',toc(timer), ...
    'actualOriginalNSEvaluations',0,'fineJacobianMatrixAllocated',false,'physicalMeshChanged',false,'newProductionSolvedAreaM2',0);
end
