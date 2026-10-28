function [update,report]=solveRocketTwoLevelPreconditioner(graph,coarse,cache,rhs,postSweeps)
%SOLVEROCKETTWOLEVELPRECONDITIONER 图扫描后作全域Galerkin粗层校正，保持固定线性映射。
% 模块一：先取得图扫描近似，再限制原完整选中分支A=-J的剩余右端。
% 模块二：固定LU解C=RAP并延拓，粗层不替换细网格，不执行非线性更新。
% 模块三：图内部残差、粗层代数残差和原NS线性精度分别记录。
if nargin<5,postSweeps=0;end
validateattributes(postSweeps,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
assert(isequaln(coarse.baseQ,cache.Q)&&isequaln(coarse.space.geometry,cache.geometry),'RocketCoarse:Cache','粗矩阵与当前完整方向缓存不匹配。');
[smoothed,graphReport]=solveRocketShellGraphPreconditioner(graph,rhs);fineAction=-applyCutCellSelectedDirection(cache,smoothed);
remaining=coarse.space.restriction*(rhs-fineAction);[correction,coarseReport]=solveRocketGalerkinCoarseSystem(coarse,remaining);
update=smoothed+coarse.space.prolongation*correction;
% 后平滑次数固定，不能依照每个右端残差自适应停机而仍称固定GMRES预条件。
for sweep=1:postSweeps
    remaining=rhs+applyCutCellSelectedDirection(cache,update);
    post=solveRocketShellGraphPreconditioner(graph,remaining);update=update+post;
end
assert(all(isfinite(update),'all'),'RocketCoarse:Nonfinite','粗层校正非有限，不裁结果。');
report=graphReport;report.status='applied-fixed-graph-scan-and-volume-Galerkin-coarse-correction';
report.coarseAlgebraicRelativeResidual=coarseReport.coarseAlgebraicRelativeResidual;
report.coarseUnknowns=coarseReport.coarseUnknowns;report.coarseGroupCount=coarseReport.coarseGroupCount;
report.coarseCorrectionModelApplications=1;report.coarseConstructionModelApplications=coarse.constructionModelApplications;
report.preconditionerResidualIsInitialGraphScan=true;report.fixedLinearMap=true;
if postSweeps>0
    report.fixedPostSweeps=postSweeps;report.coarseCorrectionModelApplications=1+postSweeps;
    report.coarseAlgebraicResidualIsBeforePostsmoothing=true;
end
end
