function [left,right,audit]=cutCellExactLineFaceLinearization(QL,QR,normal,distance,gas)
%CUTCELLEXACTLINEFACELINEARIZATION 法向线采用精确当前波系Euler导数及原冻结两点黏性项。
% 模块一：原Rusanov面导数和原近似黏性/功/导热项先按旧表达式计算。
% 模块二：仅替换Euler部分；墙、侧面与外边界的旧近似不在本函数内改变。
% 模块三：这是预条件矩阵；正式NS仍使用原完整空间算子验收。
[left,right,legacy]=cutCellLineFaceLinearization(QL,QR,normal,distance,gas);
[AL,~,~]=cutCellConservativeDerivatives(QL,normal,gas);[AR,~,~]=cutCellConservativeDerivatives(QR,normal,gas);
[EL,ER,exact]=cutCellExactEulerBranchJacobian(QL,QR,normal,gas);
AL=.5*AL;AR=.5*AR;N=size(QL,1);
for a=1:5
    AL(a,a,:)=AL(a,a,:)+reshape(.5*legacy.waveSpeedFrozen,1,1,N);
    AR(a,a,:)=AR(a,a,:)-reshape(.5*legacy.waveSpeedFrozen,1,1,N);
end
left=left-AL+EL;right=right-AR+ER;
audit=struct('kind','selected-exact-Euler-branch-plus-legacy-frozen-two-point-viscous-normal-line', ...
    'sampleRegion',exact.sampleRegion,'zeroContactVelocityCount',exact.zeroContactVelocityCount, ...
    'normalDistanceM',distance,'waveSpeedFrozen',legacy.waveSpeedFrozen,'viscousWorkVelocityFrozen',legacy.viscousWorkVelocityFrozen, ...
    'globalSmoothJacobianClaimed',false,'fullOriginalNSJacobian',false,'physicalDiscretizationReplaced',false);
end
