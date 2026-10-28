function report=testRocketMinimumStateDirection(folder)
%TESTROCKETMINIMUMSTATEDIRECTION 独立KKT证书、单位表示不变性及失败守门。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);rng(84);[U,~]=qr(randn(12,8),0);[V,~]=qr(randn(8));
A=U*diag(logspace(0,5,8))*V';truth=(1:8)';b=A*truth;
[z,r]=solveRocketMinimumStateCoefficients(A,b,.0999);assert(r.feasible);
stationary=A'*(A*z-b)+r.multiplier*z;
kkt=norm(stationary)/(norm(A,'fro')*norm(A*z-b)+r.multiplier*norm(z));
assert(kkt<1e-12&&norm(z)<norm(truth)&&r.modelDefect<=.0999+1e-10);
[scaled,rs]=solveRocketMinimumStateCoefficients(A*1e-7,b*1e-7,.0999);
unitError=norm(scaled-z)/norm(z);assert(unitError<1e-11&&rs.feasible);
[~,impossible]=solveRocketMinimumStateCoefficients([1;0],[0;1],.0999);assert(~impossible.feasible);
gas=struct('gamma',1.4,'R',287.05287);P=[.018 -1400 20 500 227;.05 -1300 -40 130 700];
Q=[P(:,1),P(:,1).*P(:,2:4),P(:,1).*(gas.R/(gas.gamma-1)*P(:,5)+.5*sum(P(:,2:4).^2,2))];
scale=[.018 28 28 28 24000];metric=buildRocketPrimitiveUpdateMetric(Q,scale,gas);u=randn(2,5);
inverse=invertRocketPrimitiveUpdateMetric(metric,applyRocketPrimitiveUpdateMetric(metric,u));inverseError=norm(inverse-u,'fro')/norm(u,'fro');assert(inverseError<1e-13);
% 只允许归约顺序舍入；任一局部残差/通量修改必须被拒绝。
rate=randn(100,5);vol=ones(100,1);net=-sum(rate,1);info=struct('netBoundaryFlux',net,'faceFlux',randn(4,5), ...
    'volumeIntegratedResidual',sum(rate,1),'balanceResidual',zeros(1,5));old=info;old.volumeIntegratedResidual=old.volumeIntegratedResidual+1e-13;
replay=assertRocketNSReplay(rate,info,rate,old,vol);bad=info;bad.faceFlux(1)=bad.faceFlux(1)+eps;
rejected=false;try,assertRocketNSReplay(rate,bad,rate,old,vol);catch e,rejected=strcmp(e.identifier,'PeakSystem:Replay');end
assert(rejected);
report=struct('status','passed-KKT-units-inverse-infeasibility-and-replay-local-rejection', ...
    'kktRelativeError',kkt,'unitRepresentationRelativeError',unitError,'primitiveMetricInverseRelativeError',inverseError, ...
    'globalReductionReplay',replay,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'minimum_state_unit.mat'),'report','A','b','z','r','-v7');
fid=fopen(fullfile(folder,'unit.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
