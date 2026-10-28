function coarse=buildRocketRefinedCoarseOperator(cache,space,parent,maximumImageBytes)
%BUILDROCKETREFINEDCOARSEOPERATOR 局部细分粗像，未细分父列原样复用，所有新子列直接重算。
% 模块一：不通过相减构造末子列，避免舍入噪声扩散成非物理稀疏非零项。
% 模块二：新像仍为同一完整细层A*P，经列归一化后固定小Gram/Cholesky。
% 模块三：构造模型作用和复制列数分账，原NS与物理网格均不改变。
if nargin<4,maximumImageBytes=256*1024^2;end
validateattributes(maximumImageBytes,{'numeric'},{'real','finite','scalar','positive'});
assert(isequaln(cache.Q,parent.baseQ)&&isequaln(cache.geometry,space.geometry)&&isequaln(parent.space.geometry,cache.geometry)&& ...
    isequaln(parent.gas,cache.gas)&&isequaln(parent.boundary,cache.boundary)&&isequaln(parent.scale,cache.scale)&&isequaln(parent.sourceIdentity,rocketSelectedDirectionSourceIdentity), ...
    'TargetCoarse:Parent','父粗像不属于当前基态、几何或计算源码。');
n=size(cache.Q,1);ng=space.groupCount;oldNg=parent.space.groupCount;T=space.restriction*parent.space.prolongation;
assert(norm(space.prolongation*T-parent.space.prolongation,Inf)<1e-12,'TargetCoarse:Nested','新空间不是父空间细分。');
parentOfGroup=zeros(ng,1);for g=1:ng,ids=find(T(g,:));assert(numel(ids)==1);parentOfGroup(g)=ids;end
childCount=accumarray(parentOfGroup,1,[oldNg 1]);d=5*ng;columns=cell(1,d);liftScale=zeros(d,1);nnzCount=0;calls=0;copied=0;timer=tic;
for b=1:5,for g=1:ng
    col=(b-1)*ng+g;pg=parentOfGroup(g);
    if childCount(pg)==1
        oldCol=(b-1)*oldNg+pg;columns{col}=parent.image(:,oldCol);liftScale(col)=parent.columnLiftScale(oldCol);copied=copied+1;
    else
        direction=zeros(n,5);vscale=1/sqrt(space.groupVolume(g));direction(space.cellToGroup==g,b)=vscale;
        value=-applyCutCellSelectedDirection(cache,direction);calls=calls+1;magnitude=norm(value(:));assert(isfinite(magnitude)&&magnitude>0,'TargetCoarse:NullImage','新子列像无效。');
        columns{col}=sparse(value(:)/magnitude);liftScale(col)=vscale/magnitude;
    end
    nnzCount=nnzCount+nnz(columns{col});assert(16*nnzCount+8*(d+1)<=maximumImageBytes,'TargetCoarse:ImageBudget','局部粗像超过原预算，不截非零。');
end,end
image=horzcat(columns{:});clear columns;buffers=whos('image');assert(buffers.bytes<=maximumImageBytes,'TargetCoarse:ImageBudget','实际稀疏像超过预算。');
G=full(image'*image);G=.5*(G+G');condition=rcond(G);assert(isfinite(condition)&&condition>1e-12,'TargetCoarse:Condition','局部粗像Gram病态，不补对角。');
[upper,flag]=chol(G);assert(flag==0,'TargetCoarse:Factor','局部粗像Gram不是数值正定。');
coarse=struct('kind','volume-based-fine-L2-minimum-residual-coarse-correction','space',space,'baseQ',cache.Q, ...
    'gas',cache.gas,'boundary',cache.boundary,'scale',cache.scale,'sourceIdentity',rocketSelectedDirectionSourceIdentity, ...
    'image',image,'columnLiftScale',liftScale,'normalMatrix',G,'upper',upper,'imageBytes',buffers.bytes, ...
    'maximumImageBytes',maximumImageBytes,'imageNonzeros',nnz(image),'reciprocalConditionEstimate',condition, ...
    'constructionModelApplications',calls,'copiedParentImageColumns',copied,'buildSeconds',toc(timer),'actualOriginalNSEvaluations',0, ...
    'fineJacobianMatrixAllocated',false,'newProductionSolvedAreaM2',0);
end
