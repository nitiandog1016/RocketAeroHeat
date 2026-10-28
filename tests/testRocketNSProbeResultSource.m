function report=testRocketNSProbeResultSource(folder)
%TESTROCKETNSPROBERESULTSOURCE 验证真实result末态来源、原root兼容和位置混用拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
source=fullfile(root,'outputs','v83_validation','riemannreservoirunit02','riemann_reservoir_unit_evidence.mat');
s=load(source,'Q','rate','result','cut','gas','boundary');a=readRocketNSProbeBase(source,'root');b=readRocketNSProbeBase(source,'result');
assert(isequaln(a.Q,s.Q)&&isequaln(a.rate,s.rate)&&isequaln(b.Q,s.result.Q)&&isequaln(b.rate,s.result.rate)&&~isequaln(a.Q,b.Q));
oldProbe=fullfile(root,'outputs','v83_validation','directnsprobeunit02','probe','direct_ns_direction_probe.mat');
c=struct('file',oldProbe,'variant','ns-probe','scale',1);[old,oldReport]=loadRocketSteadyDirectionGuess(c,s.Q,s.cut,s.gas,s.boundary,s.result.residualScale);
assert(strcmp(oldReport.sourceStateLocation,'root'));oldEvidence=load(oldProbe,'directionU');assert(isequal(old,oldEvidence.directionU));
reference=solveRocketShellSteadyKrylov(b.Q,s.cut,s.gas,s.boundary,s.result.options);
assert(reference.linearLog{1}.fullSystemTargetMet&&reference.acceptedNonlinearIterations==1);
probeFolder=fullfile(folder,'probe');p=createRocketNSDirectionProbe(probeFolder,source,source,'cut',reference.linearLog{1}.steadyDirectionEvidence.directionU,'result');
assert(p.originalNSDirectionTargetMet);c.file=fullfile(probeFolder,'direct_ns_direction_probe.mat');
[direction,loaded]=loadRocketSteadyDirectionGuess(c,b.Q,s.cut,s.gas,s.boundary,b.scale);assert(strcmp(loaded.sourceStateLocation,'result'));
options=s.result.options;options.steadyInitialGuess=c;result=solveRocketShellSteadyKrylov(b.Q,s.cut,s.gas,s.boundary,options);
assert(result.linearLog{1}.gmresSkipped&&isequaln(result.Q,reference.Q)&&isequaln(result.rate,reference.rate)&&isequaln(result.info,reference.info));
original=load(c.file);rejected=false(3,1);
for j=1:3
    bad=original;
    if j==1,bad.sourceBinding.stateLocation='root';elseif j==2,bad.sourceBinding.stateLocation='automatic';
    else,bad.sourceBinding.file=oldProbe;bad.sourceBinding.sha256=rocketStreamingFileSHA256(oldProbe);end
    f=fullfile(folder,sprintf('bad_location_%d.mat',j));save(f,'-struct','bad','-v7');config=c;config.file=f;
    try,loadRocketSteadyDirectionGuess(config,b.Q,s.cut,s.gas,s.boundary,b.scale);catch ex,rejected(j)=startsWith(ex.identifier,'SteadyGuess:');end
    assert(rejected(j));
end
report=struct('status','passed-real-result-endpoint-root-backward-compatibility-explicit-location-and-original-NS-update', ...
    'rootAndResultStatesActuallyDiffer',true,'oldRootProbeExactlyLoaded',true,'sameAcceptedStateAsFreshDefault',true, ...
    'locationRejections',rejected,'freshReferenceNSCalls',reference.residualEvaluations,'separateProbeNSCalls',p.newActualNSEvaluations, ...
    'freshMainNSCalls',result.residualEvaluations,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
Q=b.Q;gas=s.gas;boundary=s.boundary;cut=s.cut;
save(fullfile(folder,'result_source_probe_unit_evidence.mat'),'report','result','reference','Q','gas','boundary','cut','direction','loaded','source','-v7');fprintf('%s\n',jsonencode(report));
end
