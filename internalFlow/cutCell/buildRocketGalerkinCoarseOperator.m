function coarse=buildRocketGalerkinCoarseOperator(cache,space)
%BUILDROCKETGALERKINCOARSEOPERATOR 通过完整选中分支方向构造C=R*(-J)*P。
% 模块一：只保存小粗矩阵，不装配百万未知量细矩阵，不保存全部细方向列。
% 模块二：未知排列沿用MATLAB的组数x5列展开，五守恒变量缩放保持。
% 模块三：粗层无伪质量项；LU仅用于预条件校正，不能授予物理稳态资格。
assert(isequaln(space.geometry,cache.geometry),'RocketCoarse:Cache','粗空间与方向缓存的几何不一致。');
ng=space.groupCount;n=size(cache.Q,1);dimension=5*ng;C=zeros(dimension);R=space.restriction;map=space.cellToGroup;
t=tic;
for b=1:5
    for g=1:ng
        direction=zeros(n,5);direction(map==g,b)=1;action=-applyCutCellSelectedDirection(cache,direction);restricted=R*action;
        C(:,(b-1)*ng+g)=restricted(:);
    end
end
condition=rcond(C);assert(isfinite(condition)&&condition>1e-14,'RocketCoarse:Singular','粗层矩阵近奇异（rcond=%.9g），不能补对角冒充原投影。',condition);
[L,U,permutation]=lu(C,'vector');assert(all(isfinite(L),'all')&&all(isfinite(U),'all'),'RocketCoarse:Factor','粗层LU非有限。');
factorError=norm(L*U-C(permutation,:),'fro')/max(1,norm(C,'fro'));
assert(factorError<1e-12,'RocketCoarse:Factor','粗层分解不符合原矩阵。');
coarse=struct('kind','volume-Galerkin-selected-complete-spatial-Jacobian-coarse-operator','space',space,'baseQ',cache.Q, ...
    'matrix',C,'lower',L,'upper',U,'rowPermutation',permutation,'reciprocalConditionEstimate',condition, ...
    'factorizationRelativeDifference',factorError,'buildSeconds',toc(t),'constructionModelApplications',dimension, ...
    'actualOriginalNSEvaluations',0,'physicalMassOrTimeStepAdded',false,'newProductionSolvedAreaM2',0);
end
