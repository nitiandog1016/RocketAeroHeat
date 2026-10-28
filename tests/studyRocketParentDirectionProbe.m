function report=studyRocketParentDirectionProbe(folder,parentFile)
%STUDYROCKETPARENTDIRECTIONPROBE 将未合格C方向仅作为候选，重新检查完整选中J与原NS总方向。
% 模块一：C右端必须精确对应同Q最佳方向的当前剩余，不能把两个分母混用。
% 模块二：先重算完整选中J；仅模型达到原0.1时才调用原NS总方向差分。
% 模块三：原C1e-3资格不晋升，原NS/下降/温域/空间规则不变；本模块不接受新Q。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);p=load(parentFile,'rhs','value','image','detail','sourceBinding');
sourceFile=p.sourceBinding.sourceFile;assert(strcmp(rocketStreamingFileSHA256(sourceFile),p.sourceBinding.sourceFileSHA256));
s=load(sourceFile,'result','gas','boundary');r=s.result;Q=r.Q;scale=r.residualScale;baseRate=r.initialRate;rhs=baseRate./scale;
assert(r.acceptedNonlinearIterations==0&&isequaln(Q,r.initialQ));direction0=r.linearLog{end}.steadyDirectionEvidence.directionU;
savedOriginalNSDefect=r.linearLog{end}.fullSystemRelativeDefect;coarseFile=r.options.coarseOperatorFile;gas=s.gas;boundary=s.boundary;
target=r.options.linearRelativeTolerance;assert(target==.1,'ParentDirection:Target','原NS方向目标必须保持0.1。');
fdOptions=struct('differenceStep',r.options.differenceStep,'scheme',r.options.jacobianScheme,'saveEvidence',true,'spatialConservationTolerance',r.options.spatialConservationTolerance);
baseInfo=struct;for key={'primitive','gradient','viscousFaceFlux'},baseInfo.(key{1})=r.initialInfo.(key{1});end;clear r s;
t=load(coarseFile,'coarse');assert(isequaln(t.coarse.baseQ,Q)&&isequaln(t.coarse.scale,scale)&&isequaln(t.coarse.gas,gas)&&isequaln(t.coarse.boundary,boundary));cut=t.coarse.space.geometry;clear t;
sourceIdentity=rocketSelectedDirectionSourceIdentity;cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);clear baseInfo;
initialModel=applyCutCellSelectedDirection(cache,direction0);remaining=rhs+initialModel;
rhsIdentity=norm(remaining-p.rhs,'fro')/max(1,norm(p.rhs,'fro'));assert(rhsIdentity<1e-8,'ParentDirection:RHS','完整C右端不对应当前同Q增量剩余。');
delta=p.value;deltaModel=applyCutCellSelectedDirection(cache,delta);rawTotal=direction0+delta;rawModel=applyCutCellSelectedDirection(cache,rawTotal);modelCalls=3;
linearity=norm(rawModel-initialModel-deltaModel,'fro')/max(1,norm(rawModel,'fro'));assert(linearity<1e-8);
rawDefect=norm(rhs+rawModel,'fro')/norm(rhs,'fro');alpha=-dot(remaining(:),deltaModel(:))/dot(deltaModel(:),deltaModel(:));
if rawDefect<=target,totalDirection=rawTotal;totalModel=rawModel;chosenAlpha=1;
else,totalDirection=direction0+alpha*delta;totalModel=applyCutCellSelectedDirection(cache,totalDirection);modelCalls=modelCalls+1;chosenAlpha=alpha;end
chosenDefect=norm(rhs+totalModel,'fro')/norm(rhs,'fro');nsCalls=0;nsChecked=false;nsMet=false;nsDefect=NaN;slope=NaN;nsEvidence=struct;guardIdentifier='';guardMessage='';freshBaseDifference=NaN;
if chosenDefect<=target
    clear cache;
    nsCalls=1;[freshBase,~]=cutCellNavierStokesResidual(Q,cut,gas,boundary);freshBaseDifference=norm(freshBase-baseRate,'fro')/max(1,norm(baseRate,'fro'));
    assert(freshBaseDifference<1e-12,'ParentDirection:Base','原NS基态残差未复现，不执行方向检查。');
    nsCalls=nsCalls+1+double(strcmp(fdOptions.scheme,'central'));
    try
        [jv,nsEvidence]=rocketShellResidualDirectionalDerivative(Q,totalDirection,freshBase,cut,gas,boundary,scale,fdOptions);nsEvidence.value=jv;
        nsDefect=norm(rhs+jv,'fro')/norm(rhs,'fro');slope=dot(rhs(:),jv(:))/dot(rhs(:),rhs(:));upper=-1+nsDefect;
        assert(slope<=upper+512*eps*max([1 abs(slope) abs(upper)]),'ParentDirection:Descent','原NS总方向未通过原Cauchy下降界。');
        nsChecked=true;nsMet=nsDefect<=target&&slope<0;
    catch ex
        if ~(startsWith(ex.identifier,'CutNS:')||startsWith(ex.identifier,'CutExactEuler:')||startsWith(ex.identifier,'ShellJvp:')||startsWith(ex.identifier,'ParentDirection:')),rethrow(ex);end
        guardIdentifier=ex.identifier;guardMessage=ex.message;
    end
end
report=struct('status','tested-parent-C-candidate-current-selected-J-and-conditional-original-NS-total-direction-probe', ...
    'sourceBinding',p.sourceBinding,'sourceIdentity',sourceIdentity,'parentFile',parentFile,'parentFileSHA256',rocketStreamingFileSHA256(parentFile), ...
    'coarseFileSHA256',rocketStreamingFileSHA256(coarseFile),'parentCRelativeResidual',p.detail.relativeCResidual,'parentCInternalTargetMet',p.detail.internalCTargetMet, ...
    'sourceSavedOriginalNSDirectionDefect',savedOriginalNSDefect,'currentInitialModelDefect',norm(remaining,'fro')/norm(rhs,'fro'), ...
    'sameRHSIdentityDifference',rhsIdentity,'rawTotalCurrentModelDefect',rawDefect,'optimalScalar',alpha,'chosenScalar',chosenAlpha,'chosenTotalCurrentModelDefect',chosenDefect, ...
    'selectedJCombinationLinearityDifference',linearity,'parentCToSelectedNegativeJDeltaDifference',norm(p.image+deltaModel,'fro')/max(1,norm(deltaModel,'fro')), ...
    'originalNSDirectionChecked',nsChecked,'originalNSTotalDirectionDefect',nsDefect,'originalNSDirectionTarget',target,'originalNSDirectionTargetMet',nsMet, ...
    'originalNSNormalizedMeritSlope',slope,'freshOriginalNSBaseDifference',freshBaseDifference,'guardIdentifier',guardIdentifier,'guardMessage',guardMessage, ...
    'newSelectedJApplications',modelCalls,'newActualNSEvaluations',nsCalls,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalStateChanged',false,'physicalSolution',false, ...
    'originalNSResidualDenominator','original-base-rate-divided-by-original-scale','engineProcessID',feature('getpid'));
parentImage=p.image;parentRHS=p.rhs;parentCDetail=p.detail;sourceBinding=p.sourceBinding;clear p;
save(fullfile(folder,'parent_direction_probe_evidence.mat'),'Q','scale','rhs','baseRate','direction0','initialModel','remaining','delta','deltaModel','rawTotal','rawModel', ...
    'totalDirection','totalModel','nsEvidence','report','parentImage','parentRHS','parentCDetail','parentFile','sourceFile','coarseFile','sourceBinding','sourceIdentity','gas','boundary','-v7');
fid=fopen(fullfile(folder,'parent_direction_probe_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
