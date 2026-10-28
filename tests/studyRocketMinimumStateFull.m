function report=studyRocketMinimumStateFull(folder,preparationFile,candidateFile,previousFolder,segmentCount,resume)
%STUDYROCKETMINIMUMSTATEFULL 同Q24、同方向预算的原NS连续三次有界对照。
% 原方向与15维最小状态选择共享首个当前系统；每次接受后必须重建当前系统。
% 每态预先固定相同LSQR段数/15完整模型列/17候选，不依据结果临时追加。
if nargin<5,segmentCount=2;end
if nargin<6,resume=false;end
assert(ismember(segmentCount,[2 3]),'MinimumState:Budget','仅允许预声明的两段或三段预算。');
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir')||resume,'MinimumState:Overwrite','已有目录仅允许显式检查点接续。');
assert(~isfile(fullfile(folder,'full_comparison.json')),'MinimumState:Terminal','已有终态，不能重跑。');
if ~exist(folder,'dir'),mkdir(folder);end
firstP=preparationFile;firstC=candidateFile;firstPrevious=previousFolder;
branches=cell(1,2);
for method=1:2
    name={'baseline','minimum'};branchFolder=fullfile(folder,name{method});if ~exist(branchFolder,'dir'),mkdir(branchFolder);end
    steps={};timer=tic;priorElapsed=0;cachedElapsed=0;
    preparationFile=firstP;candidateFile=firstC;previousFolder=firstPrevious;
    journal=fullfile(branchFolder,'progress.json');
    if resume&&isfile(journal)
        prior=jsondecode(fileread(journal));
        if isstruct(prior)&&isscalar(prior)&&isfield(prior,'steps')
            prior.steps=asCells(prior.steps);branches{method}=prior;continue;
        end
        steps=asCells(prior);priorElapsed=sum(cellfun(@(x)x.elapsedSeconds,steps));
        assert(~isempty(steps)&&isfield(steps{end},'resultFolder')&&steps{end}.main.acceptedNonlinearIterations==1, ...
            'MinimumState:Resume','缺少已接受的完整接续检查点。');
        previousFolder=steps{end}.resultFolder;
        assert(isfile(fullfile(previousFolder,'wall_velocity_krylov_evidence.mat')));
    end
    for iteration=numel(steps)+1:3
        iterationTimer=tic;
        if iteration>1
            systemFolder=fullfile(branchFolder,sprintf('system%d',iteration));stateFile=fullfile(previousFolder,'wall_velocity_krylov_evidence.mat');
            if resume&&isfile(fullfile(systemFolder,'peak_system.json'))
                reused=loadRocketPeakSystem(fullfile(systemFolder,'peak_system.json'),stateFile);cachedElapsed=cachedElapsed+reused.constructionSeconds;
                fprintf('MINIMUM_STATE_RESUME current-system iteration=%d\n',iteration);
            else
                systemFolder=unusedFolder(systemFolder);
                prepareRocketCurrentPeakSystem(systemFolder,stateFile,fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat'),'domain');
            end
            preparationFile=fullfile(systemFolder,'peak_system.json');candidateFile='';
            for segment=1:segmentCount
                nextFolder=fullfile(branchFolder,sprintf('direction%d_%d',iteration,segment));
                if resume&&isfile(fullfile(nextFolder,'peak_zero_summary.json'))&&isfile(fullfile(nextFolder,'peak_zero_candidate_evidence.mat'))
                    saved=load(fullfile(nextFolder,'peak_zero_candidate_evidence.mat'),'preparationFileSHA256','segments','detail');
                    assert(strcmp(saved.preparationFileSHA256,rocketStreamingFileSHA256(preparationFile))&& ...
                        numel(saved.segments)==segment,'MinimumState:Resume','保存方向不属于当前系统/段次。');
                    cachedElapsed=cachedElapsed+saved.detail.elapsedSeconds;
                    fprintf('MINIMUM_STATE_RESUME completed-segment iteration=%d segment=%d\n',iteration,segment);
                else
                    nextFolder=unusedFolder(nextFolder);solveRocketPreparedPeakDirection(nextFolder,preparationFile,struct,candidateFile);
                end
                candidateFile=fullfile(nextFolder,'peak_zero_candidate_evidence.mat');
            end
        end
        p=loadRocketPeakSystem(preparationFile,fullfile(previousFolder,'wall_velocity_krylov_evidence.mat'));
        s=load(candidateFile,'value','detail','sourceBinding','segments');assert(strcmp(s.sourceBinding.sourceFileSHA256,p.sourceBinding.sourceFileSHA256));
        assert(sum([s.segments.nativeIterations])==128*segmentCount,'MinimumState:Budget','实际线性工作与预声明预算不符。');
        c=load(p.artifacts.cache.file,'cache','rate');
        [selected,selection]=selectRocketMinimumStateDirection(c.cache,s.value,c.rate,.0999);
        % 两组均计入选择计算，确保模型列和方向准备预算相同。
        direction=s.value;if method==2,direction=selected;end
        scale=c.cache.scale;rhs=c.rate./scale;clear c;
        a=load(p.artifacts.packed.file,'packed');Cvalue=applyRocketPackedGraph(a.packed,direction);clear a;
        internal=norm(rhs-Cvalue,'fro')/norm(rhs,'fro');
        stepFolder=fullfile(branchFolder,sprintf('step%d',iteration));mkdir(stepFolder);
        save(fullfile(stepFolder,'selected_direction.mat'),'direction','selected','selection','internal','candidateFile','preparationFile','-v7');
        row=struct('iteration',iteration,'method',name{method},'sourceFile',p.sourceBinding.sourceFile, ...
            'internalCRelativeResidual',internal,'internalCTarget',1e-3,'internalCTargetMet',internal<=1e-3, ...
            'selection',rmfield(selection,{'selectedModel','baselineModel','rhs'}),'probe',[],'main',[], ...
            'elapsedSeconds',0,'newProductionSolvedAreaM2',0);
        if method==2&&~selection.feasible
            row.stopReason='bounded-selection-span-infeasible';row.elapsedSeconds=toc(iterationTimer);steps{end+1}=row;break;
        end
        probeFolder=fullfile(stepFolder,'probe');
        probe=createRocketNSDirectionProbe(probeFolder,p.sourceBinding.sourceFile,p.geometryBinding.file,'domain',direction,'result');row.probe=probe;
        if ~probe.originalNSDirectionTargetMet
            row.stopReason='original-NS-direction-gate';row.elapsedSeconds=toc(iterationTimer);steps{end+1}=row;break;
        end
        resultFolder=fullfile(branchFolder,sprintf('krylov%02d',iteration));
        options=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
            'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
            'maximumStepRetries',0,'maximumLineSearch',16,'steadyInitialGuess',struct('file',fullfile(probeFolder,'direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
        row.main=continueRocketWallVelocityKrylov(resultFolder,previousFolder,options);
        row.elapsedSeconds=toc(iterationTimer);row.stopReason=row.main.stopReason;row.resultFolder=resultFolder;steps{end+1}=row;
        writeJSON(fullfile(branchFolder,'progress.json'),steps);
        if row.main.acceptedNonlinearIterations==0,break;end
        previousFolder=resultFolder;
    end
    branches{method}=struct('method',name{method},'steps',{steps},'elapsedSeconds',toc(timer)+priorElapsed+cachedElapsed, ...
        'resumedInterruptedRun',resume,'unrecordedInterruptedWorkSeconds',NaN, ...
        'elapsedSecondsKind','measured-completed-slices; interrupted-unsaved-work-excluded');
    writeJSON(fullfile(branchFolder,'progress.json'),branches{method});
end
report=struct('status','bounded-real-Q24-original-NS-comparison','branches',{branches}, ...
    'sameInitialState',firstPrevious,'firstCurrentSystemShared',true,'maximumNonlinearSteps',3, ...
    'nativeLSQRBudgetPerState',128*segmentCount,'completeModelSelectionColumnsPerState',15,'candidateBudgetPerState',17, ...
    'resumedInterruptedRun',resume,'unrecordedInterruptedWorkBudgetKnown',~resume,'speedupClaimed',false, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
writeJSON(fullfile(folder,'full_comparison.json'),report);save(fullfile(folder,'full_comparison.mat'),'report','-v7');
end
function out=asCells(value)
if iscell(value),out=value;else,out=arrayfun(@(x)x,value,'UniformOutput',false);end
end
function path=unusedFolder(base)
path=base;serial=0;
while exist(path,'dir'),serial=serial+1;path=[base '_resume' num2str(serial)];end
end
function writeJSON(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(value));fclose(fid);
end
