function [value,report]=applyRocketSpatialReduction(reduction,coarseFactor,rhs)
%APPLYROCKETSPATIALREDUCTION 固定块回代：先局部F凝聚，再粗解，最后完整F校正。
% 模块一：b_C-E*M_F*b_F是凝聚右端，x_F=M_F*(b_F-B*x_C)。
% 模块二：M_F及粗局部映射均固定；跨组F-F误差留给完整外层C检验。
assert(strcmp(reduction.kind,'five-field-spatial-local-F-approximate-Schur-reduction')&& ...
    isequal(size(rhs),[reduction.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all')&& ...
    coarseFactor.cellCount==reduction.split.coarseCells,'SpatialReduction:Input','固定空间回代或全部五量右端不符。');
split=reduction.split;bF=rhs(split.fineOriginalCVs,:);bC=rhs(split.coarseOriginalCVs,:);
first=solveRocketSpatialFBlocks(reduction.Ffactor,bF);coarseRhs=bC-applyRocketRectangularCoupling(reduction.E,first);
xC=solveRocketCoupledPatches(coarseFactor,coarseRhs);xF=solveRocketSpatialFBlocks(reduction.Ffactor,bF-applyRocketRectangularCoupling(reduction.B,xC));
value=zeros(size(rhs));value(split.fineOriginalCVs,:)=xF;value(split.coarseOriginalCVs,:)=xC;
assert(all(isfinite(value),'all'),'SpatialReduction:Nonfinite','完整空间回代非有限。');
report=struct('FLocalApplications',2,'coarseLocalApplications',1,'rectangularCouplingApplications',2,'fixedLinearMap',true, ...
    'fullFFOrCoarseInverseClaimed',false,'originalNSAcceptanceClaimed',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
