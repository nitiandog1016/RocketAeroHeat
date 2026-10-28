function report=studyRocketShellLateralCoupling(folder,sourceCase,sourceKind)
%STUDYROCKETSHELLLATERALCOUPLING 用已保存的原NS方向检查法向线遗漏的侧向耦合。
% 模块一：读取已实算的直接稳态方向及原完整NS差分，不重复进行NS调用或提交更新。
% 模块二：先补齐原冻结Rusanov侧面非对角项，再单独替换侧面Euler选中波系导数。
% 模块三：比较同一方向上的误差、五方程分量和成本；不把近似矩阵当成物理解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'LateralStudy:Overwrite','不得覆盖已有研究证据。');
assert(ismember(sourceKind,{'unit','rocket'})&&~isempty(regexp(sourceCase,'^[a-zA-Z0-9_]+$','once')));
val=fullfile(root,'outputs','v83_validation');
if strcmp(sourceKind,'unit')
    sourceFile=fullfile(val,sourceCase,'steady_newton_unit_evidence.mat');
    s=load(sourceFile,'result','cut','gas');cut=s.cut;
else
    sourceFile=fullfile(val,sourceCase,'wall_velocity_krylov_evidence.mat');
    s=load(sourceFile,'result','gas');d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;clear d;
end
r=s.result;gas=s.gas;clear s;Q=r.initialQ;scale=r.residualScale;
assert(strcmp(r.linearEquation,'steady-newton')&&numel(r.linearLog)==1&&~isempty(r.firstOperator));
e=r.linearLog{1}.steadyDirectionEvidence;assert(e.baseAcceptedIteration==0,'LateralStudy:Base','当前比较仅接受首基态对应的已保存方向。');
directionU=e.directionU;target=e.value;rhs=r.initialRate./scale;operator=r.firstOperator;
assert(isequal(size(Q),size(directionU))&&isequal(size(Q),size(target))&&all(isfinite(target),'all'));
tau=operator.pseudoTimeStepSeconds;v=cut.cells.volume(:);n=numel(v);f=cut.openFaces;
sideIDs=find(f.cellPlus>0&f.kind==3);m=double(f.cellMinus(sideIDs));p=double(f.cellPlus(sideIDs));
assert(all(floor((m-1)/operator.sourceFaces)==floor((p-1)/operator.sourceFaces)));

%% 模块四：保存线算子的实际作用；I/tau仅从预条件矩阵中代数消去。
t=tic;solved=permute(directionU,[2 3 1]);product=pagemtimes(operator.diagonal,solved);nf=operator.sourceFaces;
for layer=1:operator.layers
    ids=(layer-1)*nf+(1:nf);
    if layer>1,product(:,:,ids)=product(:,:,ids)+pagemtimes(operator.lower(:,:,ids),solved(:,:,ids-nf));end
    if layer<operator.layers,product(:,:,ids)=product(:,:,ids)+pagemtimes(operator.upper(:,:,ids),solved(:,:,ids+nf));end
end
lineAction=directionU/tau-reshape(product,5,n)';lineSeconds=toc(t);clear solved product;

%% 模块五：批量保留侧面左右邻居作用，避免装配百万未知量稀疏矩阵。
sideOffDiagonal=zeros(n,5);exactSideCorrection=zeros(n,5);exactSideFluxDerivative=zeros(numel(sideIDs),5);
batchSize=2048;rusanovSeconds=0;exactSeconds=0;largestBatchBytes=0;
for first=1:batchSize:numel(sideIDs)
    pick=first:min(first+batchSize-1,numel(sideIDs));ids=sideIDs(pick);mi=m(pick);pi=p(pick);
    normal=f.normal(ids,:);distance=sum((cut.cells.centroid(pi,:)-cut.cells.centroid(mi,:)).*normal,2);
    ql=Q(mi,:);qr=Q(pi,:);dl=directionU(mi,:).*scale;dr=directionU(pi,:).*scale;
    t=tic;[DL,DR]=cutCellLineFaceLinearization(ql,qr,normal,distance,gas);
    left=apply(DL,dl);right=apply(DR,dr);rusanovFlux=(left+right)./scale;
    % 线矩阵已含本格对角；这里只增加另一侧未知量的贡献。
    crossLeft=apply(DR,dr)./scale;crossRight=apply(DL,dl)./scale;
    for a=1:5
        sideOffDiagonal(:,a)=sideOffDiagonal(:,a)+accumarray(mi,-f.area(ids)./v(mi).*crossLeft(:,a),[n 1])+ ...
            accumarray(pi,f.area(ids)./v(pi).*crossRight(:,a),[n 1]);
    end
    rusanovSeconds=rusanovSeconds+toc(t);
    t=tic;[EL,ER]=cutCellExactLineFaceLinearization(ql,qr,normal,distance,gas);
    exactFlux=(apply(EL,dl)+apply(ER,dr))./scale;exactSideFluxDerivative(pick,:)=exactFlux;
    change=exactFlux-rusanovFlux;
    for a=1:5
        exactSideCorrection(:,a)=exactSideCorrection(:,a)+accumarray(mi,-f.area(ids)./v(mi).*change(:,a),[n 1])+ ...
            accumarray(pi,f.area(ids)./v(pi).*change(:,a),[n 1]);
    end
    exactSeconds=exactSeconds+toc(t);buffers=whos('DL','DR','EL','ER');largestBatchBytes=max(largestBatchBytes,sum([buffers.bytes]));
end
rusanovGraphAction=lineAction+sideOffDiagonal;exactSideGraphAction=rusanovGraphAction+exactSideCorrection;
comparison=struct('line',metrics(lineAction,target,rhs),'rusanovGraph',metrics(rusanovGraphAction,target,rhs), ...
    'exactSideGraph',metrics(exactSideGraphAction,target,rhs));
% 完整侧面导数替换成对散度，五量体积积分应相消；单独非对角项不具有此性质。
conservation=max(abs(sum(v.*exactSideCorrection,1))./max(1,sum(abs(v.*exactSideCorrection),1)));
assert(conservation<1e-11,'LateralStudy:Conservation','侧面完整导数替换没有保持共享通量符号。');
[~,largest]=max(max(abs(exactSideFluxDerivative),[],2));
sampleSideIndices=unique([1;round(linspace(1,numel(sideIDs),64))';largest]);
report=struct('status','measured-saved-full-NS-direction-line-versus-lateral-graph-actions', ...
    'sourceCase',sourceCase,'sourceKind',sourceKind,'cellCount',n,'sideFaceCount',numel(sideIDs), ...
    'comparison',comparison,'lineApplySeconds',lineSeconds,'rusanovSideSeconds',rusanovSeconds, ...
    'exactSideSeconds',exactSeconds,'maximumDerivativeBatchBytes',largestBatchBytes,'batchSize',batchSize, ...
    'sideReplacementVolumeConservationRelative',conservation,'sampleSideFaces',numel(sampleSideIndices), ...
    'sparseGraphStructuralNonzerosUpperBound',25*(n+2*nnz(f.cellPlus>0)), ...
    'sparseGraphCSCBytesUpperBound',16*25*(n+2*nnz(f.cellPlus>0))+8*(5*n+1), ...
    'factorizationFillMemoryMeasured',false,'newOriginalNSEvaluations',0,'newAcceptedIterations',0, ...
    'fullOriginalNSJacobian',false,'physicalEquationChanged',false,'physicalSolution',false, ...
    'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'lateral_coupling_evidence.mat'),'sourceFile','sourceCase','sourceKind', ...
    'report','lineAction','sideOffDiagonal','exactSideCorrection','rusanovGraphAction','exactSideGraphAction', ...
    'exactSideFluxDerivative','sampleSideIndices','comparison','-v7');
fid=fopen(fullfile(folder,'lateral_coupling_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function value=apply(A,d)
value=reshape(pagemtimes(A,permute(d,[2 3 1])),5,[])';
end
function out=metrics(value,target,rhs)
error=value-target;out=struct('relativeActionError',norm(error(:))/norm(target(:)), ...
    'equationRelativeActionErrors',vecnorm(error,2,1)./vecnorm(target,2,1), ...
    'steadyLinearDefect',norm(-value(:)-rhs(:))/norm(rhs(:)), ...
    'cosineToOriginalAction',dot(value(:),target(:))/(norm(value(:))*norm(target(:))));
assert(all(isfinite([out.relativeActionError out.equationRelativeActionErrors out.steadyLinearDefect out.cosineToOriginalAction])));
end
