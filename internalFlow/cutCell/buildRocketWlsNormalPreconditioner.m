function operator=buildRocketWlsNormalPreconditioner(graph,cache,maximumAdditionalBytes)
%BUILDROCKETWLSNORMALPRECONDITIONER 显式研究候选：图侧向块加完整WLS同线五带块。
% 模块一：只修正预条件矩阵，不改变原完整NS、原CV、粗像或物理边界。
% 模块二：五带法向线逐层精确块消去，原面之间仍采用原着色图近似。
% 模块三：额外持久块存储设硬预算；该数不是MATLAB进程峰值内存。
validateattributes(maximumAdditionalBytes,{'numeric'},{'real','finite','scalar','positive'});
nf=graph.sourceFaces;nl=graph.layers;n=nf*nl;extraBytes=5*25*n*8;
assert(extraBytes<=maximumAdditionalBytes,'WlsNormal:Budget','新增持久五量块超预算。');
assert(strcmp(graph.internalViscousLinearization,'frozen-two-point')&&strcmp(graph.sideEulerLinearization,'exact-branch')&& ...
    strcmp(graph.normalEulerLinearization,'exact-branch')&&isequal(graph.residualScale,cache.scale),'WlsNormal:Graph','需要当前显式精确内部Euler的原图算子。');
[delta,assembly]=assembleRocketWlsNormalCorrection(cache,maximumAdditionalBytes);operator=graph;
operator.diagonal=operator.diagonal-reshape(delta(:,:,3,:),5,5,n);
operator.lower=operator.lower-reshape(delta(:,:,2,:),5,5,n);operator.upper=operator.upper-reshape(delta(:,:,4,:),5,5,n);
operator.secondLower=-reshape(delta(:,:,1,:),5,5,n);operator.secondUpper=-reshape(delta(:,:,5,:),5,5,n);clear delta;
operator=rmfield(operator,{'lineInverse','lineUpperFactor'});
inverse=zeros(5,5,n);upper1=inverse;upper2=inverse;lower1=inverse;identity=repmat(eye(5),1,1,nf);maximumDefect=0;timer=tic;

%% 模块四：同一原面线内五带消元；不补对角、不删弱耦合或裁状态。
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);D=operator.diagonal(:,:,ids);L=operator.lower(:,:,ids);U=operator.upper(:,:,ids);
    if layer>2
        L=L-pagemtimes(operator.secondLower(:,:,ids),upper1(:,:,ids-2*nf));
        D=D-pagemtimes(operator.secondLower(:,:,ids),upper2(:,:,ids-2*nf));
    end
    if layer>1
        D=D-pagemtimes(L,upper1(:,:,ids-nf));U=U-pagemtimes(L,upper2(:,:,ids-nf));
    end
    solved=pagemldivide(D,cat(2,identity,U,operator.secondUpper(:,:,ids)));
    assert(all(isfinite(solved),'all'),'WlsNormal:Factor','五带块消去非有限，不自动补对角。');
    inverse(:,:,ids)=solved(:,1:5,:);upper1(:,:,ids)=solved(:,6:10,:);upper2(:,:,ids)=solved(:,11:15,:);lower1(:,:,ids)=L;
    maximumDefect=max(maximumDefect,max(abs(pagemtimes(D,inverse(:,:,ids))-identity),[],'all'));
end
assert(maximumDefect<1e-8,'WlsNormal:FactorDefect','五带分解没有达到原明确数值门槛。');
operator.wlsNormalFactor=struct('inverse',inverse,'upper1',upper1,'upper2',upper2,'effectiveLower1',lower1);
operator.kind='exact-selected-WLS-viscous-normal-five-band-and-existing-colored-side-graph-preconditioner';
operator.wlsNormalReport=assembly;operator.wlsNormalReport.maximumFactorDefect=maximumDefect;
operator.wlsNormalReport.factorSeconds=toc(timer);operator.wlsNormalReport.additionalPersistentBlockUpperBoundBytes=extraBytes;
operator.wlsNormalReport.maximumAdditionalBytes=maximumAdditionalBytes;
operator.fixedLinearPreconditioner=true;operator.fullOriginalNSJacobian=false;
end
