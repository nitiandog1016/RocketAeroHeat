function report=testRocketSchurFailureCheckpoint(folder)
%TESTROCKETSCHURFAILURECHECKPOINT 原数值完整回放、成功检查点及拒绝行结构化证据核对。
% 模块一：两种旧排序的所有非计时输出逐值保持，显式检查点成功路径亦相同。
% 模块二：合成近奇异2量块仍由原门槛拒绝，保存真实消元账而不返回可用因子。
% 模块三：失败部分因子被求解入口拒绝；所有测试无新J/NS或物理面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SchurObserveTest:Overwrite','不得覆盖Schur检查点测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','fieldschurunit01','field_schur_unit_evidence.mat');s=load(sourceFile,'cases');
sourceBinding=struct('sourceFile',sourceFile,'sourceFileSHA256',rocketStreamingFileSHA256(sourceFile),'physicalStateChanged',false,'unitOnly',true);
for j=1:2
    packed=s.cases{j}.packed;plain=buildRocketFieldSchurPreconditioner(packed);
    assert(isequaln(withoutTiming(plain),withoutTiming(s.cases{j}.pre)),'SchurObserveTest:Default','缺省Schur完整非计时输出未逐值复现。');
    observation=struct('folder',fullfile(folder,sprintf('assembly%d',j)),'sourceBinding',sourceBinding);
    explicit=buildRocketFieldSchurPreconditioner(packed,observation);
    assert(isequaln(withoutTiming(plain),withoutTiming(explicit)),'SchurObserveTest:Explicit','显式观察器更改了数值输出。');
    assert(exist(fullfile(observation.folder,'density_energy_failure.mat'),'file')~=2&&exist(fullfile(observation.folder,'momentum_failure.mat'),'file')~=2);
    a=load(fullfile(observation.folder,'schur_assembly_evidence.mat'));assert(isequaln(a.sourceBinding,sourceBinding));
    assert(isequaln(a.momentumMatrix.values,packed.values(2:4,2:4,:))&&isequaln(a.couplingB,packed.values(2:4,[1 5],:))&&isequaln(a.couplingE,packed.values([1 5],2:4,:)));
end
% 合成完整三行2量图，使第三个消元后对角接近奇异，前两步均是正常主元。
d=2;n=3;values=zeros(d,d,9);columns=uint32(repmat((1:3)',3,1));ptr=uint32([1;4;7;10]);diagonal=uint32([1;5;9]);
values(:,:,1)=2*eye(d);values(:,:,3)=.3*eye(d);values(:,:,5)=3*eye(d);values(:,:,6)=.2*eye(d);
values(:,:,7)=diag([.2 .15]);values(:,:,8)=.1*eye(d);target=[1 .2;.3 .06+1e-14];
values(:,:,9)=target+values(:,:,7)/values(:,:,1)*values(:,:,3)+values(:,:,8)/values(:,:,5)*values(:,:,6);
raw=struct('kind','unfactored-conservative-field-block-CSR','dimension',d,'cellCount',n,'rowPointer',ptr,'column',columns,'diagonalPosition',diagonal,'values',values);
defaultIdentifier='';try,factorRocketFieldILU0(raw);catch ex,defaultIdentifier=ex.identifier;end
assert(strcmp(defaultIdentifier,'FieldILU:FactorDefect'),'SchurObserveTest:Synthetic','合成矩阵未触发预期原逆门槛，不得改阈值。');
assemblyFile=fullfile(folder,'synthetic_assembly.mat');save(assemblyFile,'raw','sourceBinding','-v7');
observation=struct('failureFile',fullfile(folder,'synthetic_failure.mat'),'sourceBinding',sourceBinding,'assemblyFile',assemblyFile, ...
    'cellOrder',[3;1;2],'sourceFaces',1,'fieldGroup','synthetic-density-energy');
identifier='';try,factorRocketFieldILU0(raw,observation);catch ex,identifier=ex.identifier;end;assert(strcmp(identifier,defaultIdentifier));
f=load(observation.failureFile,'failure');failure=f.failure;
assert(failure.algebraicRow==3&&failure.originalControlVolume==2&&failure.originalFace==1&&failure.layer==2&& ...
    failure.maximumInverseDefect>=1e-8&&failure.originalInverseDefectThreshold==1e-8&&failure.ledgerExactlyReproducesPartialRow&&numel(failure.retainedEliminationLedger)==2&&~failure.preconditionerUsable);
rejected=false;try,solveRocketFieldILU0(failure.partialFactor,ones(3,2));catch ex,rejected=strcmp(ex.identifier,'FieldILU:Kind');end;assert(rejected);
report=struct('status','passed-default-nontiming-complete-replay-success-checkpoints-and-unchanged-original-failure-observation', ...
    'oldDefaultAndExplicitOutputsExactlyReplayed',true,'originalInverseGateUnchanged',true,'partialFactorSolveRejected',rejected, ...
    'syntheticFailureDefect',failure.maximumInverseDefect,'syntheticFailureRow',failure.algebraicRow, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'schur_checkpoint_unit_evidence.mat'),'report','sourceFile','sourceBinding','raw','failure','-v7');
fid=fopen(fullfile(folder,'schur_checkpoint_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=withoutTiming(out)
out.constructionReport=rmfield(out.constructionReport,'constructionSeconds');out.momentumFactor.factorizationReport=rmfield(out.momentumFactor.factorizationReport,'factorSeconds');out.schurFactor.factorizationReport=rmfield(out.schurFactor.factorizationReport,'factorSeconds');
end
