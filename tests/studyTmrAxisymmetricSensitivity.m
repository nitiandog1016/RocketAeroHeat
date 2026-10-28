function report=studyTmrAxisymmetricSensitivity(kind,folder,maxSteps,parentEvidenceFile)
%STUDYTMRAXISYMMETRICSENSITIVITY 官方嵌套网格上的完整五量SA稳态对照。
% 模块一：仅以冻结v66的逻辑父格复制作为数值初猜，明确不宣称保守重映射。
% 模块二：原气体、给定边界、墙热策略及固定尺度不改，逐格真实解新网格方程。
% 模块三：保存完整初末/五量账/实际成本，未稳态不得生成合格研究云图。
% 模块四：这是数值敏感性，不是实验验收、整箭求解或新增生产面积。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,kind='normal';end
if nargin<2,folder=fullfile(root,'outputs','v67_validation',char(kind));end
if nargin<3,maxSteps=90;end
if nargin<4,parentEvidenceFile='';end
validateattributes(maxSteps,{'numeric'},{'real','finite','scalar','integer','positive'});
if strcmp(kind,'outlet-profile')
    caseDefinition=createTmrAxisymmetricOutletCase(parentEvidenceFile);
else
    caseDefinition=createTmrAxisymmetricSensitivityCase(kind,parentEvidenceFile);
end
mesh=caseDefinition.mesh;gas=caseDefinition.gas;boundary=caseDefinition.boundary;
Q=caseDefinition.Q;scale=caseDefinition.residualScale;
% 数值初猜无需沿用旧稳态的BE累计账；新格点和原格点没有进行物理时间连接。
[initialRate,initialInfo]=axisymmetricRANSResidual(Q,mesh,gas,boundary);
opts=struct('maxSteps',maxSteps,'initialPseudoCFL',.5,'maximumPseudoCFL',1e6,'pseudoCFLGrowth',2, ...
    'maxNewtonIterations',24,'jacobianRefreshEvery',1,'differenceScheme','central','differenceStep',2e-6, ...
    'residualScale',scale,'steadyResidualTolerance',2e-6,'steadyRequiredSteps',3,'minimumSteadySteps',3, ...
    'conservationTolerance',5e-10,'nonlinearAbsoluteTolerance',1e-11,'nonlinearRelativeTolerance',1e-9);
target=fullfile(folder,'tmr_sensitivity_evidence.mat');
assert(~exist(target,'file'),'AxisFV:SensitivityOverwrite','禁止覆盖既有网格研究证据。');
if ~exist(folder,'dir'),mkdir(folder);end
clock=tic;implicit=advanceAxisymmetricSteady(Q,mesh,gas,boundary,opts);elapsed=toc(clock);
[rate,info]=axisymmetricRANSResidual(implicit.Q,mesh,gas,boundary);
assert(isequaln(info,implicit.finalInfo)&&max(abs(rate)./scale,[],'all')==implicit.residualNorm, ...
    'AxisFV:SensitivityReplay','实际新网格终态无法逐值重放。');
report=struct('status','actual-official-nested-grid-SA-study-not-experimental-acceptance', ...
    'kind',char(kind),'gridStudy',caseDefinition.gridStudy,'parentBinding',caseDefinition.parentBinding, ...
    'cellCount',mesh.cellCount,'initialOriginalResidual',max(abs(initialRate)./scale,[],'all'), ...
    'finalOriginalResidual',implicit.residualNorm,'steadyConverged',implicit.steadyConverged, ...
    'stopReason',implicit.stopReason,'acceptedSteps',implicit.acceptedSteps,'rejectedAttempts',implicit.rejectedAttempts, ...
    'solverResidualEvaluations',implicit.residualEvaluations,'diagnosticResidualEvaluations',2,'elapsedSeconds',elapsed, ...
    'fiveQuantityBalanceRelative',implicit.globalBalanceRelative, ...
    'signedWallHeatRangeWm2',[min(info.wall.heatFluxIntoWall),max(info.wall.heatFluxIntoWall)], ...
    'initializationOnlyNotConservativeRemap',true,'physicalTimeAdvanced',0,'pseudoTimeIsPhysical',false, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
if isfield(caseDefinition,'boundaryStudy'),report.boundaryStudy=caseDefinition.boundaryStudy;end
save(target,'report','caseDefinition','mesh','gas','boundary','Q','implicit','initialInfo');
fid=fopen(fullfile(folder,'tmr_sensitivity_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('TMR %s sensitivity: %d cells, %d accepted/%d rejected, original residual %.9g -> %.9g, steady=%d.\n', ...
    char(kind),mesh.cellCount,implicit.acceptedSteps,implicit.rejectedAttempts,report.initialOriginalResidual, ...
    report.finalOriginalResidual,implicit.steadyConverged);
end
