function report=testRocketValidatedDirectionProbe(folder)
%TESTROCKETVALIDATEDDIRECTIONPROBE 以真实小域NS方向证据检查新初猜入口和强来源拒绝。
% 小域记录来自旧已测原求解器的初始状态方向，不制造差分值或伪造流场接受历史。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','guessintegration01','steady_newton_unit_evidence.mat');s=load(sourceFile,'result','gas','boundary','cut');a=s.result;
Q=a.initialQ;scale=a.residualScale;baseRate=a.initialRate;rhs=baseRate./scale;gas=s.gas;boundary=s.boundary;cut=s.cut;
nsEvidence=a.linearLog{1}.steadyDirectionEvidence;totalDirection=nsEvidence.directionU;coarseFile=a.options.coarseOperatorFile;
sourceIdentity=rocketSelectedDirectionSourceIdentity;sourceBinding=struct('sourceFile',sourceFile,'sourceFileSHA256',rocketStreamingFileSHA256(sourceFile));
probeReport=struct('originalNSDirectionChecked',true,'originalNSDirectionTargetMet',true,'originalNSDirectionTarget',.1, ...
    'originalNSTotalDirectionDefect',a.linearLog{1}.fullSystemRelativeDefect,'physicalStateChanged',false,'newFlowAcceptances',0,'coarseFileSHA256',rocketStreamingFileSHA256(coarseFile));
validFile=fullfile(folder,'valid_probe.mat');writeProbe(validFile,probeReport);
config=struct('file',validFile,'variant','validated-probe','scale',1);[direction,loaded]=loadRocketSteadyDirectionGuess(config,Q,cut,gas,boundary,scale);
assert(isequal(direction,totalDirection)&&loaded.originalNSProbeArithmeticVerified&&~loaded.oldDerivativeReused);
options=a.options;options.steadyInitialGuess=config;result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(result.linearLog{1}.gmresSkipped&&result.linearLog{1}.fullSystemTargetMet&&result.acceptedNonlinearIterations==1&&result.residualEvaluations==3);
rejected=false(8,1);original=load(validFile);
for j=1:8
    b=original;
    if j==1,b.Q(1,1)=b.Q(1,1)*1.01;
    elseif j==2,b.scale(1)=b.scale(1)*1.01;
    elseif j==3,b.report.originalNSDirectionTargetMet=false;
    elseif j==4,b.totalDirection(1)=b.totalDirection(1)+1e-4;
    elseif j==5,b.nsEvidence.value(1)=b.nsEvidence.value(1)+1;
    elseif j==6,b.sourceIdentity.applyCutCellSelectedDirection=repmat('0',1,64);
    elseif j==7,b.sourceBinding.sourceFileSHA256=repmat('0',1,64);
    else,b.boundary.solidWall.temperature=b.boundary.solidWall.temperature+1;end
    file=fullfile(folder,sprintf('bad_probe_%02d.mat',j));save(file,'-struct','b','-v7');bad=config;bad.file=file;
    try,loadRocketSteadyDirectionGuess(bad,Q,cut,gas,boundary,scale);catch ex,rejected(j)=startsWith(ex.identifier,'SteadyGuess:');end
    assert(rejected(j));
end
report=struct('status','passed-real-small-NS-probe-loader-recomputed-arithmetic-eight-source-guards-and-fresh-main-line-search', ...
    'validDirectionExactlyPreserved',true,'rejected',rejected,'freshMainNSCalls',result.residualEvaluations,'acceptedSmallIterations',result.acceptedNonlinearIterations, ...
    'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
save(fullfile(folder,'validated_probe_unit_evidence.mat'),'report','result','validFile','config','loaded','sourceFile','gas','boundary','cut','-v7');fprintf('%s\n',jsonencode(report));
    function writeProbe(file,p)
        report=p;save(file,'Q','scale','baseRate','rhs','gas','boundary','nsEvidence','totalDirection','coarseFile','sourceIdentity','sourceBinding','report','-v7');
    end
end
