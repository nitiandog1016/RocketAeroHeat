function [inverse,certificate]=certifyRocketLocalInverse(A)
%CERTIFYROCKETLOCALINVERSE 为有界局部原块构造逆，并以原输入残差上界验收。
% 模块一：可逆行列尺度和行主元LU，不加对角项；失败返回明确不可用状态。
% 模块二：必要时补偿评价并至多两次精化，保留普通乘积误差与精确残差界的区别。
N=size(A,1);inverse=[];certificate=struct('usable',false,'upperBound',Inf,'ordinaryBLASDefect',Inf,'accurateResidualUsed',false,'refinementSteps',0);
assert(isa(A,'double')&&isreal(A)&&size(A,2)==N&&all(isfinite(A),'all'),'LocalInverse:Input','需要有限实局部方阵。');
rowMax=max(abs(A),[],2);if any(rowMax==0),return;end;R=1./rowMax;if any(~isfinite(R)),return;end
rowScaled=R.*A;columnMax=max(abs(rowScaled),[],1)';if any(columnMax==0),return;end;C=1./columnMax;if any(~isfinite(C)),return;end
balanced=rowScaled.*C';if any(~isfinite(balanced),'all')||any(A~=0&balanced==0,'all'),return;end
[L,U,p]=lu(balanced,'vector');if any(diag(U)==0)||any(~isfinite(diag(U))),return;end
right=R.*eye(N);inverse=C.*(U\(L\right(p,:)));if any(~isfinite(inverse),'all'),inverse=[];return;end
ordinary=max(abs(A*inverse-eye(N)),[],'all');gamma=N*eps/(1-N*eps);bound=ordinary+gamma*norm(A,inf)*max(abs(inverse),[],'all')/(1-gamma)^2;
if bound>=1e-8
    certificate.accurateResidualUsed=true;[residual,errorBound]=rocketCompensatedInverseResidual(A,inverse);bound=max(abs(residual),[],'all')+errorBound;
    for step=1:2
        if bound<1e-8,break;end
        right=R.*residual;proposal=inverse+C.*(U\(L\right(p,:)));[newResidual,newError]=rocketCompensatedInverseResidual(A,proposal);newBound=max(abs(newResidual),[],'all')+newError;
        if newBound>=bound,break;end
        inverse=proposal;residual=newResidual;bound=newBound;certificate.refinementSteps=step;
    end
    ordinary=max(abs(A*inverse-eye(N)),[],'all');
end
certificate.upperBound=bound;certificate.ordinaryBLASDefect=ordinary;certificate.usable=isfinite(bound)&&bound<1e-8;
end
