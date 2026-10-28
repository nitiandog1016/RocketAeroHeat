function [reference,audit]=axisymmetricJacobianLocalReference(Q,gas,scale,relativeStep)
%AXISYMMETRICJACOBIANLOCALREFERENCE 为保守量差分建立局部物理尺度，不修改状态或残差。
% 模块一：密度和总能用自身正量，动量用rho*max(声速,速度模)，避免零横向动量无扰动。
% 模块二：SA量W=rho*nuTilde用自身量级，另设64eps*mu的物理扰动底限以处理严格零值。
% 模块三：输出归一化参考及实际初始扰动审计；参考固定在当前Jacobian点，不参与微分。
% 这是明确可选的数值策略，不是边界条件、正性/精度证明，也不改变旧Jacobian默认。
gas=axisymmetricFVGas(gas);nq=size(Q,2);
assert(ismember(nq,[4 5]),'AxisLocalJacobian:State','只支持四/五守恒量。');
validateattributes(scale,{'numeric'},{'real','finite','positive','numel',nq});scale=reshape(double(scale),1,nq);
validateattributes(relativeStep,{'numeric'},{'real','finite','positive','scalar'});
P=axisymmetricFVState(Q(:,1:4),gas);rho=P(:,1);speed=max(P(:,6),vecnorm(P(:,2:3),2,2));
physical=[rho,rho.*speed,rho.*speed,abs(Q(:,4))];floorSA=[];
if nq==5
    assert(all(isfinite(Q(:,5))&Q(:,5)>=0),'AxisLocalJacobian:State','SA守恒量必须非负。');
    transport=axisymmetricFVTransport(P(:,4),gas);floorSA=64*eps*transport.mu;
    physical(:,5)=max(abs(Q(:,5)),floorSA/relativeStep);
end
reference=physical./scale;perturbation=relativeStep*reference;
assert(all(isfinite(reference)&reference>0,'all')&&all(isfinite(perturbation)&perturbation>0,'all'), ...
    'AxisLocalJacobian:Scale','局部扰动必须有限且严格正；不能用零列掩盖下溢。');
audit=struct('policy','local-conservative-acoustic-SA-roundoff-floor-v1', ...
    'relativeStep',relativeStep,'normalizedReference',reference,'initialNormalizedPerturbation',perturbation, ...
    'physicalReference',physical,'SAminimumPhysicalPerturbation',floorSA, ...
    'minimumNormalizedPerturbationByQuantity',min(perturbation,[],1), ...
    'maximumNormalizedPerturbationByQuantity',max(perturbation,[],1), ...
    'scaleHeldFixedDuringDifferentiation',true,'stateModified',false,'physicalSolution',false);
end
