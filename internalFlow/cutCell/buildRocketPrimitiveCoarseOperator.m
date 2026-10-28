function coarse=buildRocketPrimitiveCoarseOperator(cache,space,maximumImageBytes)
%BUILDROCKETPRIMITIVECOARSEOPERATOR 固定原始变量切向粗基，仍以原细层缩放L2范数拟合。
% 模块一：每个组的五种原始扰动先按当前局部Q提升，再调用已校核完整空间作用。
% 模块二：只保存有界稀疏像和小Gram因子，源Q/提升来源和计算源码同时绑定。
if nargin<3,maximumImageBytes=256*1024^2;end
validateattributes(maximumImageBytes,{'numeric'},{'real','finite','scalar','positive'});
assert(isequaln(cache.geometry,space.geometry),'PrimitiveCoarse:Geometry','粗空间和缓存几何不符。');
ng=space.groupCount;d=5*ng;columns=cell(1,d);liftScale=zeros(d,1);count=0;timer=tic;
for b=1:5,for g=1:ng
    col=(b-1)*ng+g;coeff=zeros(ng,5);vscale=1/sqrt(space.groupVolume(g));coeff(g,b)=vscale;
    direction=liftRocketPrimitiveCoarse(cache.Q,cache.gas,cache.scale,space,coeff);
    value=-applyCutCellSelectedDirection(cache,direction);magnitude=norm(value(:));assert(isfinite(magnitude)&&magnitude>0,'PrimitiveCoarse:NullImage','原始变量粗方向像无效。');
    columns{col}=sparse(value(:)/magnitude);liftScale(col)=vscale/magnitude;count=count+nnz(columns{col});
    assert(16*count+8*(d+1)<=maximumImageBytes,'PrimitiveCoarse:ImageBudget','原始变量粗像超过原内存预算，不截非零。');
end,end
image=horzcat(columns{:});clear columns;buffers=whos('image');assert(buffers.bytes<=maximumImageBytes,'PrimitiveCoarse:ImageBudget','实际原始变量粗像超过预算。');
G=full(image'*image);G=.5*(G+G');condition=rcond(G);assert(isfinite(condition)&&condition>1e-12,'PrimitiveCoarse:Condition','原始变量Gram病态，不补对角。');
[upper,flag]=chol(G);assert(flag==0,'PrimitiveCoarse:Factor','原始变量Gram不是数值正定。');
coarse=struct('kind','volume-primitive-tangent-fine-L2-minimum-residual-coarse-correction','space',space,'baseQ',cache.Q, ...
    'gas',cache.gas,'boundary',cache.boundary,'scale',cache.scale,'sourceIdentity',rocketSelectedDirectionSourceIdentity, ...
    'primitiveLiftSourceSHA256',rocketPrimitiveCoarseLiftIdentity,'image',image,'columnLiftScale',liftScale,'normalMatrix',G,'upper',upper, ...
    'imageBytes',buffers.bytes,'maximumImageBytes',maximumImageBytes,'imageNonzeros',nnz(image),'reciprocalConditionEstimate',condition, ...
    'constructionModelApplications',d,'buildSeconds',toc(timer),'actualOriginalNSEvaluations',0,'fineJacobianMatrixAllocated',false,'newProductionSolvedAreaM2',0);
end
