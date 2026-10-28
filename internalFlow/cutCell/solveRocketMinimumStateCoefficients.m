function [z,report]=solveRocketMinimumStateCoefficients(A,b,target)
%SOLVEROCKETMINIMUMSTATECOEFFICIENTS 原NS模型约束内选取最小原始变量幅度的方向。
% A的列必须对应状态度量下正交归一的真实方向；不截断状态或改变物理方程。
% 解 min ||z||，约束 ||A*z-b|| <= target*||b||；无可行方向时明确拒绝。
% 使用薄SVD和单调乘子方程，避免固定软惩罚随旧参考方向幅度漂移。
assert(isreal(A)&&isreal(b)&&all(isfinite(A),'all')&&all(isfinite(b))&&size(A,1)==numel(b), ...
    'MinimumState:Input','模型列/右端须有限且对应。');
assert(isscalar(target)&&target>0&&target<.1,'MinimumState:Target','目标须严格保留原0.1外部验收余量。');
b=b(:);den=norm(b);assert(den>0);[U,S,V]=svd(A,'econ');s=diag(S);c=U'*b;
positive=s>0;least=zeros(size(s));least(positive)=c(positive)./s(positive);
minimum=norm(A*(V*least)-b)/den;
report=struct('method','minimum-primitive-state-norm-under-original-model-discrepancy', ...
    'target',target,'minimumModelDefectInSpan',minimum,'feasible',minimum<=target, ...
    'basisColumns',size(A,2),'singularValues',s,'physicalEquationChanged',false,'newProductionSolvedAreaM2',0);
if ~report.feasible,z=zeros(size(A,2),1);return;end
lo=0;hi=max(s)^2;perpendicular=b-U*c;perp2=dot(perpendicular,perpendicular);
% 乘子使残差落在允许球边界；这不是逐次试流场的参数扫描。
for k=1:100
    lambda=(lo+hi)/2;remaining=lambda./(s.^2+lambda).*c;
    if sqrt(perp2+sum(remaining.^2))/den>target,hi=lambda;else,lo=lambda;end
end
filter=zeros(size(s));filter(positive)=s(positive)./(s(positive).^2+lo);
z=V*(filter.*c);report.multiplier=lo;report.modelDefect=norm(A*z-b)/den;report.stateMetricNorm=norm(z);
assert(report.modelDefect<=target+1e-10,'MinimumState:Residual','实际模型约束未满足。');
end
