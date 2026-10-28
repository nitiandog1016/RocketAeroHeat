function [x,report]=solveRocketRecycledGMRES(action,precondition,b,options)
%SOLVEROCKETRECYCLEDGMRES 用SVD选择多个保留方向的有界右预条件增强GMRES。
% 模块一：两遍正交化保持当前残差与保留像正交；物理方向同步消去保留分量。
% 模块二：小Arnoldi关系作SVD筛选，重新调用当前作用构造保留像，再QR归一化。
% 模块三：每周期重算完整模型残差；预算、正交或当前作用失败明确拒绝，不改原NS门槛。
if nargin<4,options=struct;end
% 研究开关：整方向/像同步缩放并自适应复正交；缺省输出与旧合同逐值一致。
robustArnoldi=false;
if isstruct(options)&&isscalar(options)&&isfield(options,'robustArnoldi')
    robustArnoldi=options.robustArnoldi;options=rmfield(options,'robustArnoldi');
    assert(islogical(robustArnoldi)&&isscalar(robustArnoldi),'RecycledGMRES:RobustOption','稳健研究开关必须为逻辑标量。');
end
% 显式失败观察器只保存原门槛拒绝时的证据；缺省配置/数值/输出合同保持。
failureEvidenceFile='';
if isstruct(options)&&isscalar(options)&&isfield(options,'failureEvidenceFile')
    failureEvidenceFile=options.failureEvidenceFile;options=rmfield(options,'failureEvidenceFile');
    assert(ischar(failureEvidenceFile)&&~isempty(failureEvidenceFile)&&exist(failureEvidenceFile,'file')~=2, ...
        'RecycledGMRES:FailureFile','失败证据必须使用未存在的明确文件。');
end
defaults=struct('restartLength',24,'maximumCycles',3,'retainedDirections',6,'maximumArnoldiSteps',72, ...
    'absoluteTolerance',1e-8,'maximumBasisBytes',1536*1024^2,'saveBasisEvidence',false, ...
    'breakdownRelativeTolerance',1e-13,'retentionRelativeThreshold',1e-12);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'RecycledGMRES:Options','未知增强Krylov选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
for key={'restartLength','maximumCycles','maximumArnoldiSteps'},validateattributes(options.(key{1}),{'numeric'},{'real','finite','scalar','integer','positive'});end
validateattributes(options.retainedDirections,{'numeric'},{'real','finite','scalar','integer','nonnegative','<',options.restartLength});
for key={'absoluteTolerance','maximumBasisBytes','breakdownRelativeTolerance','retentionRelativeThreshold'},validateattributes(options.(key{1}),{'numeric'},{'real','finite','scalar','positive'});end
assert(isa(action,'function_handle')&&isa(precondition,'function_handle')&&isnumeric(b)&&isreal(b)&&isvector(b)&&all(isfinite(b))&& ...
    islogical(options.saveBasisEvidence)&&isscalar(options.saveBasisEvidence),'RecycledGMRES:Input','作用、右端或证据选项非法。');
b=b(:);n=numel(b);m=min([options.restartLength options.maximumArnoldiSteps n]);k=min(options.retainedDirections,m-1);
% 此保守上界包括当前V/Z、旧/新U/C与构造临时量；不代表MATLAB进程峰值内存。
evidenceVectors=double(options.saveBasisEvidence)*options.maximumCycles*(2*m+4*k+12);
basisBound=8*n*(2*m+6*k+16+evidenceVectors)+8*8*(m+k+2)^2;
assert(basisBound<=options.maximumBasisBytes,'RecycledGMRES:Budget','Krylov基与保留构造预算不足，不裁向量或删CV。');
maximumActions=options.maximumArnoldiSteps+options.maximumCycles+(k+1)*max(0,options.maximumCycles-1);
x=zeros(n,1);r=b;U=zeros(n,0);C=zeros(n,0);normb=norm(b);actions=0;preconditions=0;steps=0;
history=normb;cycles={};flag=1;reason='bounded-cycles-or-Arnoldi-budget-exhausted';maximumOrthogonality=0;maximumProbe=0;
for cycle=1:options.maximumCycles
    if norm(r)<=options.absoluteTolerance,flag=0;reason='current-model-absolute-target-met';break;end
    if steps>=options.maximumArnoldiSteps,break;end
    count=min(m,options.maximumArnoldiSteps-steps);oldU=U;oldC=C;oldCount=size(C,2);startNorm=norm(r);
    if oldCount>0
        gamma=C'*r;x=x+U*gamma;r=r-C*gamma;
    end
    xStart=x;rStart=r;beta=norm(r);beforeActions=actions;beforePreconditions=preconditions;
    V=zeros(n,count+1);Z=zeros(n,count);H=zeros(count+1,count);used=0;breakdown=false;weights=[];
    if robustArnoldi,imageScales=ones(count,1);reorthPasses=zeros(count,1);localCorrelations=zeros(count,1);end
    if beta>options.absoluteTolerance
        V(:,1)=r/beta;
        for j=1:count
            z=callPrecondition(V(:,j));w=callAction(z);initialImageNorm=norm(w);
            if robustArnoldi
                % 全物理方向与其完整像同除一标量，保留所有CV和五个守恒量。
                % 零像保持缩放1，避免把有效零像方向除以realmin而产生溢出。
                if initialImageNorm>0,imageScales(j)=initialImageNorm;else,imageScales(j)=1;end
                z=z/imageScales(j);w=w/imageScales(j);initialImageNorm=norm(w);
                assert(all(isfinite(z))&&all(isfinite(w)),'RecycledGMRES:Scale','整方向缩放非有限，不能裁值。');
                for pass=1:8
                    coefficient=C'*w;w=w-C*coefficient;z=z-U*coefficient;
                    coefficient=V(:,1:j)'*w;H(1:j,j)=H(1:j,j)+coefficient;w=w-V(:,1:j)*coefficient;
                    magnitude=norm(w);localCorrelation=norm([C'*w;V(:,1:j)'*w])/max(realmin,magnitude);
                    if pass>=2&&(magnitude<=options.breakdownRelativeTolerance*max(1,initialImageNorm)||localCorrelation<=1e-12),break;end
                end
                reorthPasses(j)=pass;localCorrelations(j)=localCorrelation;Z(:,j)=z;
            else
            % 保留像的两遍投影同步作用于未知量方向，确保A*Z=V*H。
            for pass=1:2
                coefficient=C'*w;w=w-C*coefficient;z=z-U*coefficient;
            end
            Z(:,j)=z;
            for pass=1:2
                coefficient=V(:,1:j)'*w;H(1:j,j)=H(1:j,j)+coefficient;w=w-V(:,1:j)*coefficient;
            end
            end
            magnitude=norm(w);breakdown=magnitude<=options.breakdownRelativeTolerance*max(1,initialImageNorm);
            if ~breakdown,H(j+1,j)=magnitude;V(:,j+1)=w/magnitude;end
            small=H(1:j+1,1:j);target=[beta;zeros(j,1)];weights=leastSquares(small,target);
            estimate=norm(target-small*weights);history(end+1,1)=estimate; %#ok<AGROW>
            used=j;steps=steps+1;
            if estimate<=options.absoluteTolerance||breakdown,break;end
        end
        candidate=xStart+Z(:,1:used)*weights;
    else,candidate=xStart;end
    actual=b-callAction(candidate);actualNorm=norm(actual);
    if used>0,predicted=rStart-V(:,1:used+1)*(H(1:used+1,1:used)*weights);else,predicted=rStart;end
    relationError=norm(actual-predicted)/max(1,normb);
    orth=0;active=used+1-double(breakdown);
    if used>0
        orth=norm(V(:,1:active)'*V(:,1:active)-eye(active),'fro');
        if oldCount>0,orth=max(orth,norm(C'*V(:,1:active),'fro'));end
    end
    maximumOrthogonality=max(maximumOrthogonality,orth);
    if ~(orth<1e-8&&relationError<1e-8)&&~isempty(failureEvidenceFile)
        failure=struct('status','rejected-by-unchanged-original-Arnoldi-gate','orthogonalityError',orth, ...
            'modelRelationRelativeDifference',relationError,'orthogonalityThreshold',1e-8,'modelRelationThreshold',1e-8, ...
            'cycle',cycle,'arnoldiSteps',used,'modelApplications',actions,'preconditionerApplications',preconditions, ...
            'basisStorageUpperBoundBytes',basisBound,'maximumBasisBytes',options.maximumBasisBytes, ...
            'b',b,'xStart',xStart,'rStart',rStart,'candidate',candidate,'actualRemaining',actual,'predictedRemaining',predicted, ...
            'V',V(:,1:used+1),'Z',Z(:,1:used),'smallH',H(1:used+1,1:used),'weights',weights, ...
            'retainedU',oldU,'retainedC',oldC,'usableDirectionClaimed',false,'originalNSAcceptanceClaimed',false, ...
            'newProductionSolvedAreaM2',0);
        if robustArnoldi
            failure.robustArnoldi=true;failure.imageNormalizationScales=imageScales(1:used);
            failure.reorthogonalizationPasses=reorthPasses(1:used);failure.localCorrelationAfterReorthogonalization=localCorrelations(1:used);
        end
        save(failureEvidenceFile,'failure','-v7');
    end
    assert(orth<1e-8&&relationError<1e-8,'RecycledGMRES:Arnoldi','完整当前作用与Arnoldi关系或正交不符。');
    accepted=actualNorm<=startNorm+1e-10*max(1,startNorm);
    if accepted,x=candidate;r=actual;end
    history(end,1)=norm(r);
    entry=struct('cycle',cycle,'arnoldiSteps',used,'totalArnoldiSteps',steps,'retainedBefore',oldCount, ...
        'startResidualNorm',startNorm,'projectedResidualNorm',beta,'actualResidualNorm',actualNorm, ...
        'modelResidualNonincreaseVerified',accepted,'arnoldiRelationRelativeDifference',relationError, ...
        'basisOrthogonalityError',orth,'smallH',H(1:used+1,1:used),'retainedAfter',0, ...
        'retainedSingularValues',[],'retainedImageProbeRelativeDifference',NaN,'modelApplications',0,'preconditionerApplications',0);
    if robustArnoldi
        entry.imageNormalizationScales=imageScales(1:used);entry.reorthogonalizationPasses=reorthPasses(1:used);
        entry.localCorrelationAfterReorthogonalization=localCorrelations(1:used);
    end
    if options.saveBasisEvidence
        entry.V=V(:,1:used+1);entry.Z=Z(:,1:used);entry.retainedU=oldU;entry.retainedC=oldC;
        entry.xStart=xStart;entry.rStart=rStart;entry.weights=weights;entry.candidate=candidate;entry.actualRemaining=actual;
    end
    if ~accepted,reason='current-model-residual-increase-rejected';
    elseif actualNorm<=options.absoluteTolerance,flag=0;reason='current-model-absolute-target-met';
    elseif breakdown,reason='Arnoldi-breakdown-with-unresolved-current-model-residual';
    elseif cycle<options.maximumCycles&&steps<options.maximumArnoldiSteps&&k>0
        % [C,V]为像基；小矩阵D描述A*[U,Z]。SVD是当前研究策略，不冒称谐Ritz。
        D=zeros(oldCount+used+1,oldCount+used);D(1:oldCount,1:oldCount)=eye(oldCount);
        D(oldCount+1:end,oldCount+1:end)=H(1:used+1,1:used);
        [~,S,Y]=svd(D,'econ');singular=diag(S);eligible=find(singular>options.retentionRelativeThreshold*max(singular));
        ids=flipud(eligible(max(1,end-k+1):end));
        if ~isempty(ids)
            coefficients=Y(:,ids);rawU=oldU*coefficients(1:oldCount,:)+Z(:,1:used)*coefficients(oldCount+1:end,:);
            rawU=rawU./singular(ids)';rawC=zeros(n,numel(ids));
            for col=1:numel(ids),rawC(:,col)=callAction(rawU(:,col));end
            [C,R]=qr(rawC,0);assert(rcond(R)>1e-12,'RecycledGMRES:Retention','保留像相关，不补对角或裁物理向量。');U=rawU/R;
            probe=sin((1:numel(ids))'*.37);probeValue=callAction(U*probe);
            error=norm(probeValue-C*probe)/max(1,norm(probeValue));maximumProbe=max(maximumProbe,error);
            assert(error<1e-8&&norm(C'*C-eye(numel(ids)),'fro')<1e-10,'RecycledGMRES:Retention','保留提升/像的当前作用或正交失败。');
            entry.retainedAfter=numel(ids);entry.retainedSingularValues=singular(ids);entry.retainedImageProbeRelativeDifference=error;
            if options.saveBasisEvidence,entry.nextU=U;entry.nextC=C;entry.retentionCoefficients=coefficients;entry.retentionSingularValues=singular(ids);end
        else,U=zeros(n,0);C=zeros(n,0);end
    end
    entry.modelApplications=actions-beforeActions;entry.preconditionerApplications=preconditions-beforePreconditions;cycles{end+1}=entry; %#ok<AGROW>
    fprintf('RECYCLED_GMRES cycle=%d steps=%d residual=%.9g retained=%d models=%d\n',cycle,steps,norm(r)/max(realmin,normb),entry.retainedAfter,actions);
    if ~accepted||flag==0||breakdown,break;end
end
if norm(r)<=options.absoluteTolerance,flag=0;end
report=struct('method','SVD-selected-augmented-restarted-right-GMRES','flag',flag,'stopReason',reason,'options',options, ...
    'arnoldiSteps',steps,'modelApplications',actions,'preconditionerApplications',preconditions,'maximumModelApplications',maximumActions, ...
    'basisStorageUpperBoundBytes',basisBound,'maximumBasisBytes',options.maximumBasisBytes,'maximumBasisOrthogonalityError',maximumOrthogonality, ...
    'maximumRetainedImageProbeRelativeDifference',maximumProbe,'absoluteResidual',norm(r),'relativeResidual',norm(r)/max(realmin,normb), ...
    'residualHistory',history,'cycles',{cycles},'originalNSAcceptanceClaimed',false,'newProductionSolvedAreaM2',0);
if robustArnoldi
    report.robustArnoldi=true;report.imageNormalization='whole-physical-direction-and-complete-image-same-scalar';
    report.maximumReorthogonalizationPasses=8;report.localCorrelationTarget=1e-12;
    report.originalOrthogonalityThreshold=1e-8;report.originalModelRelationThreshold=1e-8;
end
    function out=callAction(v)
        assert(actions<maximumActions,'RecycledGMRES:WorkBudget','完整模型作用预算用尽。');out=action(v);actions=actions+1;
        assert(isnumeric(out)&&isreal(out)&&isequal(size(out),[n 1])&&all(isfinite(out)),'RecycledGMRES:Action','当前完整作用非法。');
    end
    function out=callPrecondition(v)
        assert(preconditions<options.maximumArnoldiSteps,'RecycledGMRES:WorkBudget','预条件作用预算用尽。');out=precondition(v);preconditions=preconditions+1;
        assert(isnumeric(out)&&isreal(out)&&isequal(size(out),[n 1])&&all(isfinite(out)),'RecycledGMRES:Precondition','预条件方向非法。');
    end
end
function x=leastSquares(A,b)
% 小矩阵数值秩处理仅用于代数最小二乘，不截细层像的非零或物理控制体。
[U,S,V]=svd(A,'econ');s=diag(S);ids=find(s>1e-13*max(s));
if isempty(ids),x=zeros(size(A,2),1);else,x=V(:,ids)*((U(:,ids)'*b)./s(ids));end
end
