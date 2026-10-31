function testSU2AsyncBinding(out,source)
%TESTSU2ASYNCBINDING 验证异步源锁/取消与研究快照版本守卫；不重新求收敛解。
% 模块一：新目录复制小型输入，启动真实进程；改动显式初始快照必须留下不稳定记录。
% 模块二：复用已有研究场测试过期编号、不同几何、伪造资格均被拒绝。
assert(~isfolder(out));mkdir(out);
for name={'case.cfg','mesh.su2'},copyfile(fullfile(source,name{1}),fullfile(out,name{1}));end
cfg=fileread(fullfile(out,'case.cfg'));cfg=regexprep(cfg,'(?m)^ITER=.*$','ITER= 30000');write(fullfile(out,'case.cfg'),cfg);
write(fullfile(out,'input_snapshot.mat'),'explicit lifecycle test input; not a MATLAB dataset');
job=embeddedSU2Job('start',out,struct('pythonExecutable','D:/anacondapython/python.exe','maximumSeconds',30));
start=tic;receipt=[];
while toc(start)<20
    [job,receipt]=embeddedSU2Job('poll',job);
    if ~isempty(receipt) && isfield(receipt,'processID'),break;end
    pause(.1);
end
assert(~job.terminal && isfield(receipt,'processID'));
write(fullfile(out,'input_snapshot.mat'),'deliberately modified after locking');
job=embeddedSU2Job('cancel',job);
while ~job.terminal && toc(start)<30,[job,receipt]=embeddedSU2Job('poll',job);pause(.1);end
assert(job.terminal && receipt.cancelled && job.supervisorExitCode==125);
stability=struct2cell(receipt.inputsStable);assert(nnz(~cellfun(@logical,stability))==1);
input=load(fullfile(source,'input_snapshot.mat'),'frozen');data=load(fullfile(source,'mapped_research.mat'),'result');
s=aeroHeatWorkflowState('init',[],struct('cfg',input.frozen.cfg));
s=aeroHeatWorkflowState('setModel',s,struct('rawMesh',input.frozen.mesh,'orientedMesh',input.frozen.mesh));
s=aeroHeatWorkflowState('setMesh',s,struct('mesh',input.frozen.mesh));
s=aeroHeatWorkflowState('confirmBoundary',s,struct('freestream',input.frozen.freestream));
s=aeroHeatWorkflowState('beginSolve',s);p=struct('result',data.result,'sourceRevision',s.solveInputRevision,'solveId',s.solveId);
bad=p;bad.solveId=bad.solveId+1;refuse(s,bad);
bad=p;bad.result.V=bad.result.V*1.1;refuse(s,bad);
bad=p;bad.result.productionQualified=true;refuse(s,bad);
bad=p;bad.result.receipt.cancelled=true;refuse(s,bad);
accepted=aeroHeatWorkflowState('completeResearch',s,p);
assert(isempty(accepted.result) && strcmp(accepted.status,'research-ready'));
report=struct('cancelTerminal',true,'changedInputDetected',true,'staleTokenRejected',true, ...
    'changedGeometryRejected',true,'qualifiedResearchRejected',true,'cancelledResultRejected',true,'newQualifiedArea',0);
write(fullfile(out,'binding_report.json'),jsonencode(report));
fprintf('Async source/input lock, owned cancellation, stale geometry/token and research qualification guards passed.\n');
end
function refuse(s,p)
failed=false;try,aeroHeatWorkflowState('completeResearch',s,p);catch,failed=true;end
assert(failed,'无效研究结果未被拒绝。');
end
function write(file,text)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
