function [value,report]=applyRocketSpatialSchur(reduction,rhs)
%APPLYROCKETSPATIALSCHUR 独立按D-E*M_F*B作用核对已装配近似Schur。
% 模块一：使用真实固定F缓存逆，含其明确的CV对角替代，不换成理想F逆。
% 模块二：该作用仅为近似Schur，不调用完整原NS。
nc=reduction.split.coarseCells;
assert(strcmp(reduction.kind,'five-field-spatial-local-F-approximate-Schur-reduction')&& ...
    isequal(size(rhs),[nc 5])&&all(isfinite(rhs),'all'),'SpatialReduction:Input','空间消元对象或粗右端非法。');
b=applyRocketRectangularCoupling(reduction.B,rhs);f=solveRocketSpatialFBlocks(reduction.Ffactor,b);
value=applyRocketRectangularCoupling(reduction.D,rhs)-applyRocketRectangularCoupling(reduction.E,f);
report=struct('rectangularCouplingApplications',3,'FLocalApplications',1,'fullFFInverseClaimed',false,'newActualNSEvaluations',0);
end
