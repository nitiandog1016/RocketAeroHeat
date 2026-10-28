function report=inspectRocketLookaheadInterfaceRejection(folder,studyFile)
%INSPECTROCKETLOOKAHEADINTERFACEREJECTION 用原内核递归定位已记录的Riemann界面越温。
% 模块一：只重放已保存的拒绝候选；核对全部单元温度已合法，物性温域不变。
% 模块二：二分实际内部面集合，原内核在单面上仍拒绝才记录，不放宽温域取样。
% 模块三：回到当前合法基态构造这些面的解析H，检查它们是否在当前选择集合内。
% 模块四：这是内部通量诊断，0新完整NS/方向/流场更新，不授予物理覆盖资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
study=jsondecode(fileread(studyFile));step=study.steps(1);
stateFile=fullfile(step.resultFolder,'wall_velocity_krylov_evidence.mat');
s=load(stateFile,'result','gas');r=s.result;gas=s.gas;clear s;
logs=r.lineLog;selected=find(strcmp({logs.errorIdentifier},'CutExactEuler:Temperature'),1,'last');assert(~isempty(selected));alpha=logs(selected).alpha;
assert(numel(r.trace)==1);trace=r.trace{1};Q=trace.initialQ;value=trace.linearUpdateU;scale=r.residualScale;clear r trace;
candidate=rocketShellTrialState(Q,value,scale,alpha,gas,'conservative');
cv=gas.R/(gas.gamma-1);rho=candidate(:,1);u=candidate(:,2:4)./rho;
T=(candidate(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);
assert(all(rho>0)&all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)));
g=load(fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat'),'domain');cut=g.domain;clear g;
f=cut.openFaces;ids=find(f.area>0&f.cellMinus>0&f.cellPlus>0);
queue={[1 numel(ids)]};bad=zeros(0,1);messages=cell(0,1);calls=0;sampledFaceInputs=0;
while ~isempty(queue)
    span=queue{end};queue(end)=[];localIDs=ids(span(1):span(2));calls=calls+1;sampledFaceInputs=sampledFaceInputs+numel(localIDs);
    try
        cutCellExactEulerFlux(candidate(f.cellMinus(localIDs),:),candidate(f.cellPlus(localIDs),:),f.normal(localIDs,:),gas);
    catch e
        if ~strcmp(e.identifier,'CutExactEuler:Temperature'),rethrow(e);end
        if numel(localIDs)==1
            bad(end+1,1)=localIDs;messages{end+1,1}=e.message; %#ok<AGROW>
        else
            middle=floor((span(1)+span(2))/2);queue{end+1}=[middle+1 span(2)];queue{end+1}=[span(1) middle]; %#ok<AGROW>
        end
    end
end
assert(~isempty(bad)&&issorted(bad));
baseRows=buildRocketInterfaceTemperatureRows(Q,cut,gas,scale,bad);
d=load(step.directionFile,'rows');selectedFace=ismember(bad,d.rows.faceIDs);originalFace=ismember(bad,step.directionReport.lookaheadSelection.originalBandFaceIDs);
baseT=baseRows.baseTemperatureK;derivative=baseRows.H*value(:);
owner=baseRows.ownerCV;neighbour=baseRows.neighbourCV;
rejectedInputs=struct('faceIDs',bad,'QL',candidate(owner,:),'QR',candidate(neighbour,:), ...
    'normal',f.normal(bad,:),'ownerCV',owner,'neighbourCV',neighbour,'gas',gas,'messages',{messages});
report=struct('status','observed-all-rejected-interior-faces-by-original-kernel-bisection-with-original-temperature-domain', ...
    'sourceStateFile',stateFile,'sourceStateSHA256',rocketStreamingFileSHA256(stateFile), ...
    'sourceDirectionSHA256',rocketStreamingFileSHA256(step.directionFile),'rejectedAlpha',alpha, ...
    'candidateCellTemperatureRangeK',[min(T) max(T)],'allCandidateCellTemperaturesInOriginalRange',true, ...
    'completeInteriorFaceCount',numel(ids),'rejectedInterfaceFaceIDs',bad, ...
    'baseInterfaceTemperatureK',baseT,'baseInterfaceSampleRegion',baseRows.sampleRegion, ...
    'currentDirectionalTemperatureChangeK',derivative,'firstOrderCandidateTemperatureK',baseT+alpha*derivative, ...
    'inCurrentSelectedConstraintSet',selectedFace,'inOriginal10KSet',originalFace, ...
    'originalKernelCalls',calls,'totalFaceInputsIncludingBisectionRechecks',sampledFaceInputs, ...
    'originalTemperatureRangeK',gas.temperatureRange,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'temperatureDomainRelaxed',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'interface_rejection.mat'),'report','baseRows','rejectedInputs','-v7');
fid=fopen(fullfile(folder,'interface_rejection.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));disp(report);
end
