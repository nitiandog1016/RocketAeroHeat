function report=testRocketDirectNSDirectionProbe(folder)
%TESTROCKETDIRECTNSDIRECTIONPROBE 原NS直接初猜的小域守门及旧无初猜完整路径回放。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
baseFile=fullfile(root,'outputs','v83_validation','riemannreservoirunit02','riemann_reservoir_unit_evidence.mat');s=load(baseFile,'Q','rate','result','cut','gas','boundary');
Q=s.Q;gas=s.gas;boundary=s.boundary;cut=s.cut;options=s.result.options;scale=s.result.residualScale;
reference=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);assert(isequaln(untimed(reference),untimed(s.result)),'DirectProbeTest:Replay','无初猜完整旧路径改变。');
probeFolder=fullfile(folder,'probe');p=createRocketNSDirectionProbe(probeFolder,baseFile,baseFile,'cut',s.result.linearLog{1}.steadyDirectionEvidence.directionU);
assert(p.originalNSDirectionTargetMet);file=fullfile(probeFolder,'direct_ns_direction_probe.mat');
config=struct('file',file,'variant','ns-probe','scale',1);[value,loaded]=loadRocketSteadyDirectionGuess(config,Q,cut,gas,boundary,scale);
assert(isequaln(value,s.result.linearLog{1}.steadyDirectionEvidence.directionU)&&~loaded.coarseImageRequired);
options.steadyInitialGuess=config;result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(result.linearLog{1}.gmresSkipped&&result.acceptedNonlinearIterations==1&&result.residualEvaluations==4&& ...
    isequaln(result.Q,s.result.Q)&&isequaln(result.rate,s.result.rate)&&isequaln(result.info,s.result.info));
original=load(file);rejected=false(8,1);
for j=1:8
    bad=original;
    if j==1,bad.Q(1)=bad.Q(1)*1.01;elseif j==2,bad.scale(1)=bad.scale(1)*1.01;
    elseif j==3,bad.boundary.solidWall.temperature=bad.boundary.solidWall.temperature+1;
    elseif j==4,bad.geometryBinding.sha256=repmat('0',1,64);elseif j==5,bad.directionU(1)=bad.directionU(1)+1e-4;
    elseif j==6,bad.nsEvidence.value(1)=bad.nsEvidence.value(1)+1;
    elseif j==7,bad.sourceIdentity.resolveNSRiemannReservoir3D=repmat('0',1,64);else,bad.sourceBinding.sha256=repmat('0',1,64);end
    badFile=fullfile(folder,sprintf('bad_probe_%02d.mat',j));save(badFile,'-struct','bad','-v7');c=config;c.file=badFile;
    try,loadRocketSteadyDirectionGuess(c,Q,cut,gas,boundary,scale);catch ex,rejected(j)=startsWith(ex.identifier,'SteadyGuess:');end
    assert(rejected(j));
end
report=struct('status','passed-direct-NS-probe-without-coarse-image-default-replay-and-eight-guards', ...
    'oldCompleteDefaultResultReplayed',true,'sameAcceptedStateAsDefault',true,'freshMainNSCalls',result.residualEvaluations, ...
    'separateProbeNSCalls',p.newActualNSEvaluations,'separateDefaultReplayNSCalls',reference.residualEvaluations,'rejected',rejected, ...
    'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
save(fullfile(folder,'direct_ns_probe_unit_evidence.mat'),'report','result','Q','gas','boundary','cut','config','loaded','-v7');fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});
end
