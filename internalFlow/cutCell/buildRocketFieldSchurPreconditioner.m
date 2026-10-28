function pre=buildRocketFieldSchurPreconditioner(packed,observation)
%BUILDROCKETFIELDSCHURPRECONDITIONER 构造三动量/密度能量的五量近似Schur预条件。
% 模块一：C=[A B;E D]按守恒量分组；S≈D-E*diag(A)^{-1}*B，完整五量均保留。
% 模块二：三动量A用原模式ILU，密度/能量S用完整二跳模式ILU；无时间正则或补对角。
% 模块三：只改变预条件器，不改C、原NS、边界或物理网格；先核对内存预算。
% 模块四：显式检查点在分解前保存组装矩阵；原门槛失败时另存结构化证据。
if nargin<2,observation=struct([]);end
observe=~isempty(observation);
if observe
    assert(isstruct(observation)&&isscalar(observation)&&all(ismember(fieldnames(observation),{'folder','sourceBinding'}))&& ...
        all(isfield(observation,{'folder','sourceBinding'}))&&ischar(observation.folder)&& ...
        isstruct(observation.sourceBinding)&&isscalar(observation.sourceBinding)&&~exist(observation.folder,'dir'), ...
        'FieldSchur:Observation','检查点需要新目录和明确来源，不能覆盖已有证据。');
    mkdir(observation.folder);
end
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    ~isfield(packed,'enrichmentReport')&&~isfield(packed,'flowReorderingReport'),'FieldSchur:Kind','需要原始完整五量C。');
n=packed.cellCount;nb=numel(packed.column);owners=repelem((1:n)',diff(double(packed.rowPointer)));timer=tic;
P=spones(sparse(owners,double(packed.column),1,n,n));productBound=full(sum(P,1))*full(sum(P,2));
assert(64*(productBound+2*nb)+128*n<=packed.options.maximumConstructionBytes,'FieldSchur:Budget','Schur符号数组估算超预算。');
pattern=spones(P*P);clear P;ns=nnz(pattern);assert(ns<2^32);
persistentBytes=nb*21*8+ns*4*8+(3*nb+ns)*4+3*(n+1)*4+3*n*4+n*13*8+2*n*8;
constructionBytes=nb*25*8+2*persistentBytes+(nb+ns)*48+n*256;
assert(persistentBytes<=packed.options.maximumPersistentBytes&&constructionBytes<=packed.options.maximumConstructionBytes, ...
    'FieldSchur:Budget','五量Schur子块/因子或构造数组估算超预算，不裁连接。');
fprintf('FIELD_SCHUR cells=%d base=%d schur=%d persistentBound=%.0f constructionEstimate=%.0f\n',n,nb,ns,persistentBytes,constructionBytes);
[columns,rowOwners]=find(pattern');clear pattern;
counts=accumarray(rowOwners,1,[n 1]);sptr=uint32([1;1+cumsum(counts)]);position=sparse(columns,rowOwners,(1:ns)',n,n);sdiag=uint32(full(diag(position)));clear rowOwners;
scalar=[1 5];momentum=2:4;localInverse=zeros(3,3,n);maximumLocalDefect=0;
for row=1:n
    D=packed.values(momentum,momentum,packed.diagonalPosition(row));condition=rcond(D);
    assert(isfinite(condition)&&condition>0,'FieldSchur:LocalSingular','局部动量块奇异，不补对角。');
    inv=D\eye(3);defect=max(abs(D*inv-eye(3)),[],'all');assert(all(isfinite(inv),'all')&&defect<1e-8,'FieldSchur:LocalDefect','局部动量逆未通过原1e-8。');
    localInverse(:,:,row)=inv;maximumLocalDefect=max(maximumLocalDefect,defect);
end
B=packed.values(momentum,scalar,:);E=packed.values(scalar,momentum,:);S=zeros(2,2,ns);
for first=1:4096:nb
    ids=first:min(first+4095,nb);positions=full(position(sub2ind([n n],double(packed.column(ids)),owners(ids))));assert(all(positions>0));
    S(:,:,positions)=packed.values(scalar,scalar,ids);
end
for row=1:n
    for entry=double(packed.rowPointer(row)):double(packed.rowPointer(row+1))-1
        pivot=double(packed.column(entry));ids=double(packed.rowPointer(pivot)):double(packed.rowPointer(pivot+1))-1;
        positions=full(position(sub2ind([n n],double(packed.column(ids)),repmat(row,numel(ids),1))));assert(all(positions>0));
        left=E(:,:,entry)*localInverse(:,:,pivot);S(:,:,positions)=S(:,:,positions)-pagemtimes(left,B(:,:,ids));
    end
    if mod(row,50000)==0,fprintf('FIELD_SCHUR assembly=%d/%d\n',row,n);end
end
clear position localInverse;
M=struct('kind','unfactored-conservative-field-block-CSR','dimension',3,'cellCount',n,'rowPointer',packed.rowPointer, ...
    'column',packed.column,'diagonalPosition',packed.diagonalPosition,'values',packed.values(momentum,momentum,:));
if observe
    momentumMatrix=M;schurMatrix=struct('kind','unfactored-conservative-field-block-CSR','dimension',2,'cellCount',n, ...
        'rowPointer',sptr,'column',uint32(columns),'diagonalPosition',sdiag,'values',S);
    couplingB=B;couplingE=E;algebraicOrder=packed.order;sourceBinding=observation.sourceBinding;
    assemblyReport=struct('basePatternBlocks',nb,'twoHopSchurPatternBlocks',ns,'maximumLocalMomentumInverseDefect',maximumLocalDefect, ...
        'persistentFactorStorageUpperBoundBytes',persistentBytes,'algorithmArrayConstructionEstimateBytes',constructionBytes, ...
        'checkpointOnly',true,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    sourceIdentity=struct('buildRocketFieldSchurPreconditioner',rocketStreamingFileSHA256(which('buildRocketFieldSchurPreconditioner')), ...
        'factorRocketFieldILU0',rocketStreamingFileSHA256(which('factorRocketFieldILU0')));
    assemblyFile=fullfile(observation.folder,'schur_assembly_evidence.mat');
    save(assemblyFile,'momentumMatrix','schurMatrix','couplingB','couplingE','algebraicOrder','sourceBinding','assemblyReport','sourceIdentity','-v7');
    clear momentumMatrix schurMatrix couplingB couplingE algebraicOrder;
    configuration=struct('failureFile',fullfile(observation.folder,'momentum_failure.mat'),'sourceBinding',sourceBinding, ...
        'assemblyFile',assemblyFile,'cellOrder',packed.order,'sourceFaces',packed.sourceFaces,'fieldGroup','three-momentum');
    M=factorRocketFieldILU0(M,configuration);
    save(fullfile(observation.folder,'momentum_factor_evidence.mat'),'M','sourceBinding','assemblyFile','sourceIdentity','-v7');
else,M=factorRocketFieldILU0(M);end
S=struct('kind','unfactored-conservative-field-block-CSR','dimension',2,'cellCount',n,'rowPointer',sptr,'column',uint32(columns), ...
    'diagonalPosition',sdiag,'values',S);clear columns;
if observe
    configuration.failureFile=fullfile(observation.folder,'density_energy_failure.mat');configuration.fieldGroup='density-total-energy';
    S=factorRocketFieldILU0(S,configuration);
else,S=factorRocketFieldILU0(S);end
pre=struct('kind','fixed-five-conservation-momentum-density-energy-approximate-Schur','cellCount',n,'order',packed.order, ...
    'rowPointer',packed.rowPointer,'column',packed.column,'owners',uint32(owners),'B',B,'E',E,'momentumFactor',M,'schurFactor',S, ...
    'constructionReport',struct('basePatternBlocks',nb,'twoHopSchurPatternBlocks',ns,'maximumLocalMomentumInverseDefect',maximumLocalDefect, ...
    'persistentFactorStorageUpperBoundBytes',persistentBytes,'algorithmArrayConstructionEstimateBytes',constructionBytes, ...
    'sourceOperatorAndLibraryWorkspaceExcluded',true,'processPeakMemoryClaimed',false,'constructionSeconds',toc(timer), ...
    'allOriginalCellsAndFiveEquationsKept',true,'physicalEquationChanged',false,'fixedLinearMap',true,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0));
end
