function coarse=buildRocketLeastSquaresCoarseOperator(cache,space,maximumImageBytes)
%BUILDROCKETLEASTSQUARESCOARSEOPERATOR 以完整细层残差范数构造稀疏像空间粗校正。
% 模块一：真实CV体积定义粗基，保存A*P的稀疏列，不构造细层Jacobian。
% 模块二：像列归一化后用小Gram矩阵固定Cholesky；不得补对角隐瞒秩亏。
% 模块三：R_V仍用于守恒审计；拟合采用原缩放细层L2范数，区别于Galerkin粗残差。
if nargin<3,maximumImageBytes=256*1024^2;end
validateattributes(maximumImageBytes,{'numeric'},{'real','finite','scalar','positive'});
assert(isequaln(space.geometry,cache.geometry),'RocketCoarse:Cache','细层缓存与粗空间几何不一致。');
ng=space.groupCount;n=size(cache.Q,1);dimension=5*ng;columns=cell(1,dimension);liftScale=zeros(dimension,1);count=0;t=tic;
for b=1:5,for g=1:ng
    column=(b-1)*ng+g;volumeScale=1/sqrt(space.groupVolume(g));direction=zeros(n,5);direction(space.cellToGroup==g,b)=volumeScale;
    image=-applyCutCellSelectedDirection(cache,direction);magnitude=norm(image(:));assert(isfinite(magnitude)&&magnitude>0,'RocketCoarse:NullImage','粗方向的完整空间像为零或非有限。');
    columns{column}=sparse(image(:)/magnitude);liftScale(column)=volumeScale/magnitude;count=count+nnz(columns{column});
    assert(16*count+8*(dimension+1)<=maximumImageBytes,'RocketCoarse:ImageBudget','粗像空间超过预先规定内存预算，不删除非零项。');
end,end
image=horzcat(columns{:});clear columns;buffers=whos('image');assert(buffers.bytes<=maximumImageBytes,'RocketCoarse:ImageBudget','实际稀疏像存储超过预算。');
G=full(image'*image);G=.5*(G+G');condition=rcond(G);
assert(isfinite(condition)&&condition>1e-12,'RocketCoarse:NormalCondition','归一化粗像Gram矩阵病态（rcond=%.9g），不补对角。',condition);
[upper,flag]=chol(G);assert(flag==0&&all(isfinite(upper),'all'),'RocketCoarse:NormalFactor','粗像Gram不是数值正定。');
coarse=struct('kind','volume-based-fine-L2-minimum-residual-coarse-correction','space',space,'baseQ',cache.Q, ...
    'gas',cache.gas,'boundary',cache.boundary,'scale',cache.scale,'sourceIdentity',rocketSelectedDirectionSourceIdentity, ...
    'image',image,'columnLiftScale',liftScale,'normalMatrix',G,'upper',upper,'imageBytes',buffers.bytes, ...
    'maximumImageBytes',maximumImageBytes,'imageNonzeros',nnz(image),'reciprocalConditionEstimate',condition, ...
    'constructionModelApplications',dimension,'buildSeconds',toc(t),'actualOriginalNSEvaluations',0, ...
    'fineJacobianMatrixAllocated',false,'newProductionSolvedAreaM2',0);
end
