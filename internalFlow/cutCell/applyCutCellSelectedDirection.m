function [value,report]=applyCutCellSelectedDirection(cache,directionU)
%APPLYCUTCELLSELECTEDDIRECTION 应用已校核的完整选中分支空间方向，不调用原NS。
% 模块一：缓存特征迹传播到WLS约束和开放通量，未冻结应有的边界状态敏感性。
% 模块二：完整WLS黏性/导热/功链式项与Euler相加；不提供非线性接受资格。
% 模块三：求解器须用原NS重新核对最终组合方向及实际候选，不能只凭本值验收。
Q=cache.Q;gas=cache.gas;scale=cache.scale;dq=directionU.*scale;rho=Q(:,1);u=Q(:,2:4)./rho;cv=gas.R/(gas.gamma-1);
T=(Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);du=(dq(:,2:4)-u.*dq(:,1))./rho;
dT=(dq(:,5)-sum(u.*dq(:,2:4),2)+(.5*sum(u.^2,2)-cv*T).*dq(:,1))./(rho*cv);dP=[dq(:,1) du dT];
given=zeros(numel(cache.geometry.openFaces.area),5);ids=cache.outsideFaceIDs;
given(ids,:)=reshape(pagemtimes(cache.boundaryJacobian,permute(dP(cache.outsideOwner,:),[2 3 1])),5,[])';
bd=struct('primitiveDerivative',given,'incoming',cache.state.incoming,'wall',cache.state.wall,'boundarySensitivityIncluded',true);
viscous=cutCellWlsViscousDirection(Q,directionU,cache.baseInfo,cache.geometry,gas,cache.boundary,scale,cache.gradient,cache.state,bd);
euler=cutCellSelectedEulerDirection(Q,directionU,cache.baseInfo,cache.geometry,gas,cache.boundary,scale,cache.state,bd,cache.waves);
value=euler+viscous;assert(all(isfinite(value),'all'),'CutDirection:Nonfinite','缓存完整方向作用非有限。');
report=struct('scheme','selected-analytic','actualScalarStep',NaN,'actualOriginalResidualEvaluations',0, ...
    'selectedBranchOnly',true,'boundaryStateSensitivityIncluded',true,'fullWLSViscousAndWorkChainIncluded',true, ...
    'nonlinearAcceptanceClaimed',false,'physicalTimeAdvanced',0,'newProductionSolvedAreaM2',0);
end
