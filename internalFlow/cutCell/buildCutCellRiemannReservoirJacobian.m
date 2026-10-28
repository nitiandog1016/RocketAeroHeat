function out=buildCutCellRiemannReservoirJacobian(Q,cut,gas,boundary)
%BUILDCUTCELLRIEMANNRESERVOIRJACOBIAN 批量构造精确储库采样迹对域内五原始量的导数。
% 模块一：给定储库不微分；由精确压力根、波扇及选中接触侧传播真实采样状态。
% 模块二：sample左侧守恒导数乘dQ/dP，给缓存及非缓存WLS/开放通量共同使用。
% 模块三：固定当前波系/入流选择，不声称跨接触或全局光滑；原NS仍独立验收。
state=buildCutCellRiemannReservoirBoundaryState(Q,cut,gas,boundary);ids=state.faceID;N=numel(ids);
out=struct('faceID',ids,'slot',state.slot,'primitiveJacobian',zeros(5,5,N),'basePrimitive',state.primitive,'incoming',state.incoming, ...
    'sampleRegion',state.sampleRegion,'selectedBranchOnly',true,'newOriginalNSEvaluations',0);
if N==0,return;end
f=cut.openFaces;owner=max(f.cellMinus(ids),f.cellPlus(ids));normal=f.normal(ids,:);pick=f.cellMinus(ids)==0;normal(pick,:)=-normal(pick,:);
r=state.reservoirPrimitive;cv=gas.R/(gas.gamma-1);QR=[r(:,1),r(:,1).*r(:,2:4),r(:,1).*(cv*r(:,5)+.5*sum(r(:,2:4).^2,2))];
[~,~,~,sample]=cutCellExactEulerBranchJacobian(Q(owner,:),QR,normal,gas);
assert(norm(sample.basePrimitive-state.primitive,'fro')/max(1,norm(state.primitive,'fro'))<1e-12, ...
    'CutDirection:RiemannBase','储库导数没有复现实际批量采样迹。');
rho=Q(owner,1);u=Q(owner,2:4)./rho;T=(Q(owner,5)-.5*rho.*sum(u.^2,2))./(rho*cv);K=zeros(5,5,N);
K(1,1,:)=1;
for a=1:3
    K(a+1,1,:)=reshape(u(:,a),1,1,N);K(a+1,a+1,:)=reshape(rho,1,1,N);K(5,a+1,:)=reshape(rho.*u(:,a),1,1,N);
end
K(5,1,:)=reshape(cv*T+.5*sum(u.^2,2),1,1,N);K(5,5,:)=reshape(rho*cv,1,1,N);
out.primitiveJacobian=pagemtimes(sample.left,K);assert(all(isfinite(out.primitiveJacobian),'all'));
end
