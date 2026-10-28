function report=studyRiemannReservoirSelectedDirection(folder)
%STUDYRIEMANNRESERVOIRSELECTEDDIRECTION 在真实Q20核对新储库完整解析方向及差分敏感性。
% 模块一：先精确回放原NS基态，重建只属于当前Q20的完整WLS/Euler/储库采样缓存。
% 模块二：同一旧种子方向对比解析模型和三种原NS差分步，定位失效而不继承方向资格。
% 模块三：保存真实当前缓存供新方向构造；不更新Q、不增加合格面积或拼接历史。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','krylov20','wall_velocity_krylov_evidence.mat');s=load(sourceFile,'result','gas','boundary');
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;
Q=s.result.Q;gas=s.gas;boundary=s.boundary;scale=s.result.residualScale;direction=s.result.linearLog{1}.steadyDirectionEvidence.directionU;
[rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);assert(isequaln(rate,s.result.rate)&&isequaln(info,s.result.info));clear s;
timer=tic;cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);buildSeconds=toc(timer);
[value,detail]=applyCutCellSelectedDirection(cache,direction);rhs=rate./scale;analyticDefect=norm(rhs+value,'fro')/norm(rhs,'fro');
steps=[1e-6 5e-7 2.5e-7];fdValues=cell(3,1);fdEvidence=cell(3,1);differences=zeros(3,1);defects=differences;
for j=1:3
    [fdValues{j},fdEvidence{j}]=rocketShellResidualDirectionalDerivative(Q,direction,rate,cut,gas,boundary,scale, ...
        struct('differenceStep',steps(j),'scheme','forward','saveEvidence',true,'spatialConservationTolerance',1e-9));
    differences(j)=norm(value-fdValues{j},'fro')/max(1,norm(value,'fro'));
    defects(j)=norm(rhs+fdValues{j},'fro')/norm(rhs,'fro');
end
assert(all(differences<2e-5),'RiemannDirection:Real','真实完整解析方向与原NS差分未达到既定数值一致性。');
sourceIdentity=rocketSelectedDirectionSourceIdentity;
report=struct('status','tested-current-Q20-complete-reservoir-selected-direction-and-three-original-NS-steps', ...
    'sourceFile',sourceFile,'sourceSHA256',rocketStreamingFileSHA256(sourceFile),'sourceStateLocation','result', ...
    'currentBaseQExactlyReplayed',true,'sourceFaceCount',size(cut.inputGeometry.F,1),'cellCount',size(Q,1), ...
    'analyticOriginalRHSRelativeDefect',analyticDefect,'originalNSRelativeDefects',defects,'originalNSDifferenceSteps',steps, ...
    'analyticToOriginalNSRelativeDifferences',differences,'directionOriginalNSPointOneMet',defects(1)<=.1, ...
    'cacheBuildSeconds',buildSeconds,'newActualNSEvaluations',4,'newSelectedModelApplications',1, ...
    'newFlowAcceptances',0,'physicalStateChanged',false,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
save(fullfile(folder,'riemann_selected_direction_real_evidence.mat'),'Q','gas','boundary','scale','direction','rate','info','cache','value','detail', ...
    'fdValues','fdEvidence','report','sourceIdentity','-v7');
fid=fopen(fullfile(folder,'riemann_selected_direction_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
