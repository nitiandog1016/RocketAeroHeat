function report=studyRocketWlsDirection(folder,sourceCase)
%STUDYROCKETWLSDIRECTION 在真实krylov07失败方向上核对完整WLS及选中分支空间导数。
% 模块一：释放旧图缓存，只沿保存方向作解析分解与中心/半步原NS差分。
% 模块二：分别统计两点图、完整黏性/功速度、开放边界敏感性，未知不填零。
% 模块三：这是方向一致性研究，不提交新Q、不授予稳态/物理/面积资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RocketWlsDirection:Overwrite','不得覆盖真实方向研究。');
assert(strcmp(sourceCase,'krylov07'),'RocketWlsDirection:Source','本轮明确使用已定位的krylov07失败方向。');
sourceFile=fullfile(root,'outputs','v83_validation',sourceCase,'wall_velocity_krylov_evidence.mat');
s=load(sourceFile,'result','gas','boundary','report');Q=s.result.initialQ;scale=s.result.residualScale;gas=s.gas;boundary=s.boundary;
directionU=s.result.linearLog{1}.steadyDirectionEvidence.directionU;savedForwardAction=s.result.linearLog{1}.steadyDirectionEvidence.value;
baseInfo=pick(s.result.initialInfo);sourceSolverActualCalls=s.report.totalNewEquationResidualEvaluations;clear s;
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;clear d;
state=cutCellBoundaryStencilState(Q,cut,gas,boundary);before=memory;t=tic;cache=buildCutCellGradientOperator(cut,state);cacheBuildSeconds=toc(t);after=memory;
t=tic;bd=cutCellBoundaryDirectionalState(Q,directionU,cut,gas,boundary,scale,state);boundaryDirectionSeconds=toc(t);
t=tic;[viscous,vd]=cutCellWlsViscousDirection(Q,directionU,baseInfo,cut,gas,boundary,scale,cache,state,bd);viscousSeconds=toc(t);
t=tic;[euler,ed]=cutCellSelectedEulerDirection(Q,directionU,baseInfo,cut,gas,boundary,scale,state,bd);eulerSeconds=toc(t);analytic=euler+viscous;
frozenViscous=cutCellWlsViscousDirection(Q,directionU,baseInfo,cut,gas,boundary,scale,cache,state,[]);
refinements=cell(2,1);errors=zeros(2,4);h0=1e-6/max(abs(directionU),[],'all');calls=0;
for refinement=1:2
    h=h0/2^(refinement-1);t=tic;
    [rp,full]=cutCellNavierStokesResidual(Q+h*directionU.*scale,cut,gas,boundary);ip=pick(full);clear full;
    [rm,full]=cutCellNavierStokesResidual(Q-h*directionU.*scale,cut,gas,boundary);im=pick(full);clear full;calls=calls+2;elapsed=toc(t);
    fullAction=(rp-rm)./scale/(2*h);g=(ip.gradient-im.gradient)/(2*h);
    fv=(ip.viscousFaceFlux-im.viscousFaceFlux)/(2*h);fc=(ip.convectiveFaceFlux-im.convectiveFaceFlux)/(2*h);
    errors(refinement,:)=[relative(analytic,fullAction),relative(vd.gradientDerivative,g),relative(vd.faceViscousFluxDerivative,fv),relative(ed.faceEulerFluxDerivative,fc)];
    assert(max(errors(refinement,:))<1e-3,'RocketWlsDirection:FD','真实保存方向的原NS/梯度/通量差分检查未通过，不能视为可用导数。');
    refinements{refinement}=struct('h',h,'ratePlus',rp,'rateMinus',rm,'gradientPlus',ip.gradient,'gradientMinus',im.gradient, ...
        'faceViscousPlus',ip.viscousFaceFlux,'faceViscousMinus',im.viscousFaceFlux,'faceEulerPlus',ip.convectiveFaceFlux, ...
        'faceEulerMinus',im.convectiveFaceFlux,'elapsedSeconds',elapsed);
    clear ip im;
end
prior=load(fullfile(root,'outputs','v83_validation','graphdirection01','graph_direction_evidence.mat'),'graphAction');
graphRelative=norm(prior.graphAction(:)-fullAction(:))/norm(fullAction(:));analyticRelative=norm(analytic(:)-fullAction(:))/norm(fullAction(:));
nf=size(cut.inputGeometry.F,1);nl=numel(cut.cells.volume)/nf;layerGraphError=zeros(nl,1);layerAnalyticError=layerGraphError;
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);layerGraphError(layer)=norm(prior.graphAction(ids,:)-fullAction(ids,:),'fro');
    layerAnalyticError(layer)=norm(analytic(ids,:)-fullAction(ids,:),'fro');
end
report=struct('status','passed-actual-real-failed-direction-full-WLS-and-selected-complete-NS-half-step-check', ...
    'sourceCase',sourceCase,'cellCount',size(Q,1),'fullNSMaxRelativeDifferences',errors(:,1), ...
    'gradientMaxRelativeDifferences',errors(:,2),'viscousFluxMaxRelativeDifferences',errors(:,3),'eulerFluxMaxRelativeDifferences',errors(:,4), ...
    'analyticL2RelativeDifference',analyticRelative,'graphL2RelativeDifference',graphRelative, ...
    'equationL2RelativeDifferences',vecnorm(analytic-fullAction,2,1)./vecnorm(fullAction,2,1), ...
    'savedForwardActionL2Difference',norm(savedForwardAction(:)-fullAction(:))/norm(fullAction(:)), ...
    'layerGraphErrorL2',layerGraphError,'layerAnalyticErrorL2',layerAnalyticError, ...
    'workVelocityChainNorm',norm(vd.workVelocityU(:)),'frozenBoundaryViscousDifferenceNorm',norm(viscous(:)-frozenViscous(:)), ...
    'cacheBuildSeconds',cacheBuildSeconds,'boundaryDirectionSeconds',boundaryDirectionSeconds,'viscousDirectionSeconds',viscousSeconds,'eulerDirectionSeconds',eulerSeconds, ...
    'matlabBytesBeforeCache',before.MemUsedMATLAB,'matlabBytesAfterCache',after.MemUsedMATLAB,'peakMemoryMeasured',false, ...
    'sourceSolverActualCalls',sourceSolverActualCalls,'extraValidationOriginalNSEvaluations',calls,'sourceSolverHistoryChanged',false, ...
    'newAcceptedIterations',0,'steadyConverged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
allEvidence={struct('directionU',directionU,'boundaryDirection',bd,'viscous',viscous,'viscousDetail',vd,'euler',euler,'eulerDetail',ed,'refinements',{refinements})};
mkdir(folder);save(fullfile(folder,'wls_direction_real_evidence.mat'),'report','sourceFile','Q','baseInfo','cut','gas','boundary','scale','state','allEvidence', ...
    'frozenViscous','savedForwardAction','-v7');
fid=fopen(fullfile(folder,'wls_direction_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=pick(info)
out=struct;for key={'primitive','gradient','viscousFaceFlux','convectiveFaceFlux','wallFlux','wallPressure'},out.(key{1})=info.(key{1});end
end
function error=relative(a,b)
error=max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end
