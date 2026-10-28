function coarse=buildRocketMixedCoarseOperator(spatial,primitive,remaining,options)
%BUILDROCKETMIXEDCOARSEOPERATOR 从已验证像筛选固定互补方向，不拼接超预算的完整两套像。
% 模块一：保留全部保守空间列；以保存剩余上的Schur投影收益对原始变量列排序。
% 模块二：每个新增候选重新核对完整小Gram正定性/条件数、独立性和真实存储预算。
% 模块三：筛选只在构造时执行，后续右端变化不重新选列，预条件仍为固定线性映射。
if nargin<4,options=struct;end
defaults=struct('maximumImageBytes',256*1024^2,'maximumAdditionalColumns',48, ...
    'minimumOrthogonalImageSquaredNorm',1e-8,'minimumRelativeGain',1e-10);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'MixedCoarse:Options','未知混合粗空间选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumAdditionalColumns,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
for key={'maximumImageBytes','minimumOrthogonalImageSquaredNorm','minimumRelativeGain'},validateattributes(options.(key{1}),{'numeric'},{'real','finite','scalar','positive'});end
assert(isequaln(spatial.baseQ,primitive.baseQ)&&isequaln(spatial.space.geometry,primitive.space.geometry)&& ...
    isequaln(spatial.gas,primitive.gas)&&isequaln(spatial.boundary,primitive.boundary)&&isequaln(spatial.scale,primitive.scale)&& ...
    isequaln(spatial.sourceIdentity,rocketSelectedDirectionSourceIdentity)&&isequaln(primitive.sourceIdentity,spatial.sourceIdentity)&& ...
    strcmp(primitive.primitiveLiftSourceSHA256,rocketPrimitiveCoarseLiftIdentity),'MixedCoarse:Sources','两套输入不属于同一基态、原几何、边界物性或计算/提升来源。');
assert(isequal(size(remaining),size(spatial.baseQ))&&all(isfinite(remaining),'all'),'MixedCoarse:RHS','筛选剩余非法。');
S=spatial.image;P=primitive.image;ns=size(S,2);np=size(P,2);timer=tic;
buffers=whos('S');assert(buffers.bytes<=options.maximumImageBytes,'MixedCoarse:Budget','保守底空间已超过预算。');
G0=full(S'*S);G0=.5*(G0+G0');GP=full(P'*P);GP=.5*(GP+GP');cross=full(S'*P);
assert(norm(G0-spatial.normalMatrix,'fro')<1e-10*max(1,norm(G0,'fro'))&&rcond(G0)>1e-12,'MixedCoarse:Base','保守像或小因子来源不符。');
[U0,flag]=chol(G0);assert(flag==0,'MixedCoarse:Base','底空间小Gram非正定。');
Y=U0\(U0'\cross);H=GP-cross'*Y;H=.5*(H+H');vector=remaining(:);
baseCoefficients=U0\(U0'\(S'*vector));baseRemaining=vector-S*baseCoefficients;c=P'*baseRemaining;
selected=zeros(0,1);excluded=false(np,1);reasons=repmat({''},np,1);history=zeros(0,7);nnzCount=nnz(S);G=G0;upper=U0;stop='maximum-additional-columns';
for attempt=1:np
    if numel(selected)>=options.maximumAdditionalColumns,break;end
    if isempty(selected),diagonal=diag(H);correlation=c;
    else
        small=H(selected,selected);solution=small\H(selected,:);diagonal=diag(H)-sum(H(:,selected).*solution',2);
        correlation=c-H(:,selected)*(small\c(selected));
    end
    valid=~excluded;valid(selected)=false;independent=diagonal>options.minimumOrthogonalImageSquaredNorm;
    newly=valid&~independent;excluded(newly)=true;if any(newly),reasons(newly)={'dependent-at-selection-tolerance'};end
    score=-inf(np,1);valid=valid&independent;
    for j=find(valid)'
        bytes=16*(nnzCount+nnz(P(:,j)))+8*(ns+numel(selected)+2);
        if bytes>options.maximumImageBytes,excluded(j)=true;reasons{j}='image-memory-bound';else,score(j)=correlation(j)^2/diagonal(j);end
    end
    [gain,j]=max(score);
    if ~isfinite(gain),stop='no-independent-candidate-within-memory-bound';break;end
    if gain<=options.minimumRelativeGain*max(1,sum(vector.^2)),stop='additional-projection-gain-below-construction-threshold';break;end
    trial=[selected;j];candidate=[G0 cross(:,trial);cross(:,trial)' GP(trial,trial)];candidate=.5*(candidate+candidate');condition=rcond(candidate);
    [factor,flag]=chol(candidate);
    if flag~=0||~isfinite(condition)||condition<=1e-12,excluded(j)=true;reasons{j}='full-Gram-conditioning-or-positive-definiteness';continue;end
    selected=trial;G=candidate;upper=factor;nnzCount=nnzCount+nnz(P(:,j));
    right=[S'*vector;P(:,selected)'*vector];weights=upper\(upper'\right);
    predicted=vector-S*weights(1:ns)-P(:,selected)*weights(ns+1:end);
    bytes=16*nnzCount+8*(ns+numel(selected)+1);
    history(end+1,:)=[numel(selected) j diagonal(j) gain condition bytes norm(predicted)/max(1,norm(vector))]; %#ok<AGROW>
end
if numel(selected)<options.maximumAdditionalColumns&&strcmp(stop,'maximum-additional-columns'),stop='candidate-attempt-budget-exhausted';end
image=[S P(:,selected)];buffers=whos('image');assert(buffers.bytes<=options.maximumImageBytes,'MixedCoarse:Budget','实际联合稀疏像超过预算，不截非零。');
assert(norm(full(image'*image)-G,'fro')<1e-10*max(1,norm(G,'fro')),'MixedCoarse:Gram','筛选后实际像与小Gram不符。');
primitiveSpace=rmfield(primitive.space,'geometry');
coarse=struct('kind','mixed-conservative-primitive-fine-L2-minimum-residual-coarse-correction','space',spatial.space,'primitiveSpace',primitiveSpace, ...
    'baseQ',spatial.baseQ,'gas',spatial.gas,'boundary',spatial.boundary,'scale',spatial.scale,'sourceIdentity',spatial.sourceIdentity, ...
    'liftSourceIdentity',rocketMixedCoarseLiftIdentity,'selectedPrimitiveColumns',selected,'baseConservativeColumns',ns, ...
    'columnLiftScale',[spatial.columnLiftScale;primitive.columnLiftScale(selected)],'image',image,'normalMatrix',G,'upper',upper, ...
    'imageBytes',buffers.bytes,'imageNonzeros',nnz(image),'maximumImageBytes',options.maximumImageBytes,'reciprocalConditionEstimate',rcond(G), ...
    'selectionOptions',options,'selectionRemaining',remaining,'selectionHistory',history,'excludedCandidateReasons',{reasons}, ...
    'selectionStopReason',stop,'constructionModelApplications',0,'copiedSourceImageColumns',size(image,2),'buildSeconds',toc(timer), ...
    'actualOriginalNSEvaluations',0,'fineJacobianMatrixAllocated',false,'newProductionSolvedAreaM2',0);
end
